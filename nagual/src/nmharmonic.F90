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
module nmharmonic

  use nmparameter
  use nmmath
  use nmbasis
  use nmshell

  implicit none

    public :: harmonic_real_normalized
    public :: harmonic_xyz_to_spherical
    public :: harmonic_ssxyz_integral
    public :: harmonic_pointer_size

    private

contains

  subroutine harmonic_xyz_to_spherical(x,y,z,r,t,p)
  ! Transform Cartesian to Spherical coordinates. 
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    real(8) :: p,r,t,x,y,z

    ! Radius 
    r = sqrt( x**2 + y**2 + z**2 )   

    ! Elevation angle 
    if (r.lt.ntolnum) then
      t = 0.0
    else
      t = acos(z/r)
    end if

    ! Azimuth angle 
    if (abs(x).lt.ntolnum) then
      if (abs(y).lt.ntolnum) then
        p = 0.0
      else if (y.lt.0.0) then
        p = 1.5*pi
      else
        p = 0.5*pi
      end if
    else if (x.gt.0.0) then
      p = atan(y/x)
    else 
      p = atan(y/x) + pi
    end if
  end subroutine

  subroutine harmonic_real_normalized(lmax,theta,phi,rsh)
  ! Calculate real spherical harmonics.
  ! W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery,
  ! Numerical Recipes, 2nd Edition, Cambridge University Press
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    integer :: lmax
    real(8) :: phi,theta
    real(8) :: rsh(*)

    integer :: fl,l,m
    real(8) :: norm,normzero,x
    real(8) :: p(0:nmaxlk,0:nmaxlk),s(0:nmaxlk),c(0:nmaxlk)

    ! Calculate argument for Legendre polynomials 
    x = cos(theta)

    ! 1. Special case THETA = 0 
    if (x.eq.1.0) then
      p(:,:) = 0.0
      p(:,0) = 1.0
    ! 2. Special case THETA = Pi 
    else if (x.eq.-1.0) then
      p(:,:) = 0.0
      p(0,0) = 1.0
      do l=1,lmax
        p(l,0) = -p(l-1,0)
      end do
    ! 3. General case ==> Evaluate over recurrence relations 
    else
      s(1) = sqrt((1.0-x)*(1.0+x))
      do l=2,lmax
        s(l) = s(l-1)*s(1)
      end do

      do l=0,lmax
        m = l
        if (m.eq.0) then
          p(l,m) = 1.0
        else
          p(l,m) = s(m)*math_factoriald(2*m-1)
          m = l - 1
          p(l,m) = x*(2*m+1)*p(l-1,m)
          if (l.gt.1) then
            do m=0,l-2
              p(l,m) = (x*(2*l-1)*p(l-1,m) - (l+m-1)*p(l-2,m))/(l-m)
            end do
          end if
        end if
      end do
    end if

    ! Normalization 
    do l=0,lmax
      normzero = sqrt(float(2*l+1)/(2.0*pi))
      p(l,0) = normzero*p(l,0)
      do m=1,l
        norm = sqrt(math_factorial(l-m)/math_factorial(l+m))*normzero
        p(l,m) = norm*p(l,m)
      end do
    end do

    ! Calculate trigonometric functions 
    if (lmax.gt.0) then
      if (phi.eq.0.0) then
        s(:) = 0.0
        c(:) = 1.0
      else
        s(1) = sin(phi)
        c(1) = cos(phi)
        do m=2,lmax
          s(m) = s(1)*c(m-1)+c(1)*s(m-1)
          c(m) = c(1)*c(m-1)-s(1)*s(m-1)
        end do
      end if
    end if

    ! Build the normalized Real Spherical Harmonics 
    do l=0,lmax
      fl = l*l+l+1
      rsh(fl) = p(l,0)/sqrt(2.0)
      do m=1,l
        rsh(fl-m) = p(l,m)*s(m)
        rsh(fl+m) = p(l,m)*c(m)
      end do
    end do
  end subroutine

  real(8) function harmonic_ssxyz_integral(la,ma,lb,mb,i,j,k)
  ! Calculation of angular ECP integrals.
  ! < S(LA,MA) | S(LB,MB) | x^i y^j z^k >
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    integer :: i,j,k,la,lb,ma,mb
    ! Description of some variables
    ! ctostm: Cartesian TO Spherical Transformation Matrices.
    ! ptostm: Polynom TO Spherical harmonics Transformation Matrices.
    integer al,am,ax,ay,az,bl,bm,lc
    real(8) :: norm

    ! Angularity of the polynom 
    lc = i + j + k

    ! If polynom is not present use orthonormality 
    if (lc.eq.0) then
      if ((la.eq.lb).and.(ma.eq.mb)) then
        harmonic_ssxyz_integral = 1.0
      else
        harmonic_ssxyz_integral = 0.0
      end if
    ! If polynom is present 
    else
      ! Relabel to use the smallest sum 
      if (la.ge.lb) then
        bl = la
        am = ma
        al = lb
        bm = mb
      else
        bl = lb
        am = mb
        al = la
        bm = ma
      end if
      harmonic_ssxyz_integral = 0.0
      ! Sum coefficient products if required 
      if (lc+al.ge.bl) then
        do ax=0,al
          do ay=0,al-ax
            az = al - ax - ay
            norm = sqrt(math_factoriald(2*al+1)/(4.0*pi*& 
                   math_factoriald(2*ax-1)*math_factoriald(2*ay-1)*&
                   math_factoriald(2*az-1)))
            harmonic_ssxyz_integral = harmonic_ssxyz_integral + norm*&
            ctostm(bm,raop(ax,ay,az),al)*ptostm(am,gaop(ax+i,ay+j,az+k),bl)
          end do
        end do
      end if
    end if
  end function

  integer function harmonic_pointer_size(l)
  ! Dimension for spherical harmonics pointer work fields.
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    integer :: l

    harmonic_pointer_size = (l + 1)**2
  end function

end module
