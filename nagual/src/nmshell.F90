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
module nmshell
  
  use nmparameter
  use nmtypes
  use nmfile
  use nmmath

  implicit none

    public :: shell_print
    public :: shell_copy
    public :: shell_load
    public :: shell_nco
    public :: shell_gto_radius
    public :: shell_initialize_pointers
    public :: shell_dimgaop

    integer,public :: gaop(0:nmaxli+2,0:nmaxli+2,0:nmaxli+2)
    integer,public :: raop(0:nmaxli+2,0:nmaxli+2,0:nmaxli+2)

    private

contains

  subroutine shell_load(s,n,l,k,z,c)
  ! Load basis set for a given shell
  ! Roberto Flores-Moreno, Jul 2008, Aug 2008
  implicit none
    type(nshell) :: s
    integer :: l,k,n
    real(8) :: z(*),c(*)

    integer :: i

    ! Load values
    s%n = n  ! principal quantum number for basis and radial power for ECPs
    s%l = l  ! Angularity
    s%k = k  ! Contraction order
    do i=1,k
      s%z(i) = z(i)  ! Exponents
      s%c(i) = c(i)  ! Coefficients
    end do
  end subroutine 

  subroutine shell_copy(src,dst)
  ! Copy shell from src to dst
  ! Roberto Flores-Moreno, Jul 2008-Aug 2008
  implicit none
    type(nshell) :: src,dst

    dst%n = src%n
    dst%l = src%l
    dst%k = src%k
    dst%z = src%z
    dst%c = src%c
    dst%norm = src%norm
  end subroutine 

  subroutine shell_print(tape,s)
  ! Print shell basis to tape
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape
    type(nshell) :: s

    integer :: i
    real(8) :: norm,t

    write(tape,'(t2,i2,x,i2,3x,i3)') s%n,s%l,s%k
    norm = 2.0**s%l*(2.0/pi)**0.75
    do i=1,s%k
      t = s%c(i)
      t = t/(norm*(s%z(i)**(2*s%l+3))**0.25)
      write(tape,'(t2,f30.10,x,f30.10)') s%z(i),t
    end do
  end subroutine 

  integer function shell_nco(l)
  ! Determines number of cartesian orbitals in shell
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: l

    shell_nco = ((l+1)*(l+2))/2
  end function

  integer function shell_dimgaop(l)
  ! Dimension for generalized atomic orbital pointer
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: l

    shell_dimgaop = (l**3 - l)/6 + l**2 + 2*l + 1
  end function

  subroutine shell_initialize_pointers
  ! Initialize pointers
  ! Roberto Flores-Moreno, 2008
  implicit none

    integer :: rcount,gcount,l,lx,ly,lz

    gcount = 0
    do l=0,nmaxli+2
      rcount = 0
      do lx=l,0,-1
        do ly=l-lx,0,-1
          lz = l - lx - ly
          rcount = rcount + 1
          gcount = gcount + 1
          raop(lx,ly,lz) = rcount
          gaop(lx,ly,lz) = gcount
        end do
      end do
    end do
  end subroutine

  real(8) function shell_gto_radius(z,d,n,t)
  ! Generate atomic radii for a primitive Gaussian.
  ! Roberto Flores-Moreno, 2004, 2018
  implicit none
    integer :: maxnri
    parameter (maxnri = 40)

    integer :: n
    real(8) :: d,t,z
    ! Description of some variables
    ! d: Coefficient.
    ! n: Radial power.
    ! t: Treshold for radius.
    ! z: Exponent.

    integer :: i
    real(8) :: delta,dg,g,ldt,zr

    ldt = log(abs(d)/abs(t))
    ! Analytical expression for n = 0 
    if (n.eq.0) then
      shell_gto_radius = max(ldt/z,ntolnum)
    ! Find outer radius by Newton-Raphson procedure 
    else
      shell_gto_radius = max(ldt/z,ntolnum)
      ! For large radial powers the guess should be corrected 
      if (n.gt.0) then
        g = max(sqrt(float(n)/(2*abs(d))),ntolnum)
        if (g.gt.shell_gto_radius) then
          shell_gto_radius = 0.5*(shell_gto_radius+g)
        end if
      end if
      do i=1,maxnri
        zr = z*shell_gto_radius
        g = ldt + float(n)*log(shell_gto_radius) - zr*shell_gto_radius
        dg = float(n)/shell_gto_radius - 2*zr
        delta = g/dg
        shell_gto_radius = max(shell_gto_radius-delta,ntolnum)
        if (abs(delta).lt.1.0e-10) then
          shell_gto_radius = shell_gto_radius**2
          return
        end if
      end do
      ! Reaching this point is failure
      call file_error('shell: No outer radius found in gto_radius')
    end if
  end function

end module
