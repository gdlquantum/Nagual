! This file is part of Nagual software.
!
!    Nagual is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, version 3 of the License.
!
!    Nagual is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Nagual.  If not, see <http://www.gnu.org/licenses/>.
!
! Cite as follows:
! 
! R. Flores-Moreno, H. N. Gonzalez-Ramirez, J. F. H. Lew-Yee, J. M. del Campo,
! J. J. Villalobos-Castro, J. A. Flores-Ramos, J. A. Guerrero-Cruz,
! B. Zuniga-Gutierrez, Nagual 1, Guadalajara Jal., Mexico (2020)
!
!###################################################
!   Nagual: Multicomponent many body calculations.
!   Copyright (C) 2006-2020 Nagual developers.
!
!   Contact: r.flores@academicos.udg.mx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nmvector

  use nmparameter
  use nmtypes
  use nmfile
  use nmmath

  implicit none

    public :: vector_constructor
    public :: vector_copy
    public :: vector_is_zero
    public :: vector_add
    public :: vector_substract
    public :: vector_norm
    public :: vector_normalize
    public :: vector_scale
    public :: vector_ortho
    public :: vector_colineal
    public :: vector_parallel
    public :: vector_dot
    public :: vector_distance
    public :: vector_cross
    public :: vector_rotate
    public :: vector_reflect

    private

contains

  subroutine vector_constructor(this,x,y,z)
  ! Set values of vector components
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)
    real(8) :: x,y,z

    this(1) = x
    this(2) = y
    this(3) = z
  end subroutine

  subroutine vector_copy(src,dest)
  ! Copy vector information
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: src(3),dest(3)

    dest(1:3) = src(1:3)
  end subroutine

  subroutine vector_ortho(this,other)
  ! Get other vector orthogonal to this
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),other(3)

    call vector_constructor(other,1.0,0.0,0.0)
    if (vector_is_zero(this)) return

    if (vector_colineal(this,other)) then
      call vector_constructor(other,0.0,1.0,0.0)
    else
      call vector_substract(other,this)
      call vector_normalize(other)
    end if
  end subroutine

  logical function vector_colineal(this,other)
  ! Check if these two vectors are colineal
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),other(3)

    real(8) :: dot
    real(8) :: tmp(3)

    call vector_cross(this,other,tmp)
    if (vector_is_zero(tmp)) then
      vector_colineal = .true.
    else
      vector_colineal = .false.
    end if
  end function

  logical function vector_parallel(this,other)
  ! Check if these two vectors are paralel
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),other(3)

    real(8) :: tmp(3)

    vector_parallel = .false.
    if (vector_colineal(this,other)) then
      if (vector_dot(this,other).gt.0.0) then
        vector_parallel = .true.
      end if
    end if
  end function

  real(8) function vector_dot(a,b)
  ! Evaluates scalar product of these two vectors
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: a(3),b(3)

    vector_dot = sum(a(:)*b(:))
  end function

  real(8) function vector_distance(a,b)
  ! Evaluates distance between points a and b
  ! Roberto Flores-Moreno, 2018
  implicit none
    real(8) :: a(3),b(3)

    real(8) :: x,y,z

    x = b(1) - a(1)
    y = b(2) - a(2)
    z = b(3) - a(3)
    vector_distance = sqrt(x*x + y*y + z*z)
  end function

  subroutine vector_cross(a,b,c)
  ! Evaluatex a x b cross product
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: a(3),b(3),c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine

  subroutine vector_add(this,other)
  ! Add other to this
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),other(3)

    this(1:3) = this(1:3) + other(1:3)
  end subroutine

  subroutine vector_substract(this,other)
  ! Substract other from this
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),other(3)

    this(1:3) = this(1:3) - other(1:3)
  end subroutine

  subroutine vector_normalize(this)
  ! Normalize this vector
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)

    real(8) :: norm

    if (vector_is_zero(this)) then
      call file_error('Unable to normalize vector')
    else
      norm = vector_norm(this)
      this(1:3) = this(1:3)/norm
    end if
  end subroutine

  subroutine vector_scale(this,sf)
  ! Scale elements of this vector
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)
    real(8) :: sf

    this(1:3) = this(1:3)*sf
  end subroutine

  subroutine vector_rotate(this,axis,angle)
  ! Rotate this vector around axis by angle 
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),axis(3)
    real(8) :: angle  ! in degrees

    real(8) :: rotmat(3,3)

    call build_rotation_matrix(axis,angle,rotmat)
    this = matmul(rotmat,this)
  end subroutine

  subroutine build_rotation_matrix(this,phi,rotmat)
  ! Build rotation matrix around a given axis
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)
    real(8) :: phi  ! in degrees
    real(8) :: rotmat(3,3)

    real(8) :: x,y,z,cp,sp,ct,norm

    if (.not.vector_is_zero(this)) then
      norm = vector_norm(this)
      x = this(1)/norm
      y = this(2)/norm
      z = this(3)/norm
      cp = cos(phi*pi/180.0)
      sp = sin(phi*pi/180.0)
      ct = 1.0 - cp
      rotmat(1,1) = x*x*ct + cp
      rotmat(1,2) = x*y*ct - z*sp
      rotmat(1,3) = x*z*ct + y*sp
      rotmat(2,1) = y*x*ct + z*sp
      rotmat(2,2) = y*y*ct + cp
      rotmat(2,3) = y*z*ct - x*sp
      rotmat(3,1) = z*x*ct - y*sp
      rotmat(3,2) = z*y*ct + x*sp
      rotmat(3,3) = z*z*ct + cp
    else
      rotmat(1,1) = 1.0
      rotmat(1,2) = 0.0
      rotmat(1,3) = 0.0
      rotmat(2,1) = 0.0
      rotmat(2,2) = 1.0
      rotmat(2,3) = 0.0
      rotmat(3,1) = 0.0
      rotmat(3,2) = 0.0
      rotmat(3,3) = 1.0
    end if
  end subroutine

  logical function vector_is_zero(this)
  ! Check if the vector norm is negligible
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)

    real(8) :: sp

    vector_is_zero = .false.
    if (vector_norm(this).lt.ntolnum) vector_is_zero = .true.
  end function

  real(8) function vector_norm(this)
  ! Get norm of a vector
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3)

    vector_norm  = sqrt(sum(this(1:3)**2))
  end function

  subroutine vector_reflect(this,normal)
  ! Reflect vector on a plne whose normal is given
  ! Roberto Flores-Moreno (2010)
  implicit none
    real(8) :: this(3),normal(3)

    real(8) :: s
    real(8) :: other(3)

    if (.not.vector_is_zero(normal)) then
      call vector_copy(normal,other)
      call vector_normalize(other)
      s = sum(this(1:3)*other(1:3))
      this(1:3) = this(1:3) - 2.0*s*other(1:3)
    else
      this(1:3) = 0.0
    end if
  end subroutine

end module 
