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
module nmbessel
! Literature:
!
! Half‐numerical evaluation of pseudopotential integrals
! R. Flores‐Moreno et al. Journal of computational chemistry 27 (9), 1009-1019
!
! L.E. McMurchie, E. Davidson, J. Comp. Phys. 44, 289 (1981)
!
! G. Arfken, H.J. Weber, Mathematical Methods for Physicists,
! Fourth Edition, Academic Press London, 1995, pp. 688

  use nmparameter 
  use nmfile

  implicit none

    public :: bessel_initialize
    public :: bessel_k

    private
      integer :: zabscis
      parameter (zabscis = 16*100) 

      real(8) :: ktab(0:nmaxlk+6,0:zabscis),kcoef(nmaxlk+6)  ! 6 because of Taylor series pretab

contains

  subroutine bessel_k(n,z,k)
  ! Calculate exponential weighted modified spherical Bessel 
  ! function K(z) for the semi-local ECP integrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    real(8) :: eps
    parameter (eps = 1.E-7)

    integer :: n
    real(8) :: z
    real(8) :: k(0:nmaxlk)

    integer :: i,j,l,m,loop
    real(8) :: fa,fb
    real(8) :: a(0:nmaxlk+6),d(0:nmaxlk+6)

    ! 1. z = 0 ==> Zero order approximation 
    if (z.le.eps) then
      if (z.le.0.0) then
        k(0)= 1.0 
        do i=1,n
          k(i) = 0.0 
        end do
      else
        k(0) = 1.0 - z 
        do i=1,n
          k(i) = k(i-1)*z/float(i+i+1)
        end do
      end if
    ! 2. 0 < z < 16 ==> Taylor series expansion 
    else if (z.lt.16.0) then
      loop = n + 5
      fb = float(zabscis/16)
      i = nint(z*fb)
      fa = z - float(i)/fb
      fb = 1.0
      d(0:loop) = ktab(0:loop,i)
      k(0:n) = d(0:n)
      do i=1,5
        j = loop - i 
        ! (I-1)th derivatives 
        m = j + 1
        a(0:m) = d(0:m)
        ! Calculate Ith derivatives
        d(0) = a(1) - a(0)
        do l=1,j
          m = l + 1
          d(l) = kcoef(l)*(a(l-1)-a(m))-a(l)+a(m)
        end do
        ! Add term to funcion values
        fb = fb*fa/float(i)
        k(0:n) = k(0:n) + fb*d(0:n)
      end do
    ! 3. z > 16 ==> Asymptotic formula 
    else 
      a(0) = 0.5/z
      k(0:n) = a(0) 
      do l=1,n
        fa = float((l+1)*l)
        do i=1,l-1
          k(l) = k(l) + fa*a(i)
          fa = fa*float((l+i+1)*(l-i))
        end do
        a(l) = -a(0)*a(l-1)/float(l)
        k(l) = k(l) + fa*a(l)
      end do
    end if
  end subroutine

  subroutine bessel_initialize
  ! Calculation and tabulation of K(z) for Taylor series expansion in 
  ! bessel_k. The K(z) values are calculated by a truncated series expansion.
  ! Roberto Flores-Moreno, 2004, 2018
  implicit none
    integer :: maxterm
    parameter (maxterm = 200)

    integer :: i,l,m,n
    real(8) :: u,ztab
    real(8) :: t(0:maxterm),f(0:maxterm+nmaxlk+6)

    ! Initialization 
    ktab(:,:) = 0.0
    ktab(0,0) = 1.0

    ! Loop over all points of K(z) 
    do i=1,zabscis 
      ztab = float(i)/float(zabscis/16) 
      ! Calculation of K(0) values 
      n = 0
      u = ztab*ztab/2.0
      f(n) = 1.0
      t(n) = exp(-ztab)
      l = int(0.25*sqrt(1.0+16.0*u))
      ! Series expansion 
      do while ((t(n)/f(n).gt.ntolnum).or.(n.le.l))
        ktab(0,i) = ktab(0,i) + t(n)/f(n)
        n = n + 1
        if (n.gt.maxterm) then
          call file_error('bessel_initialize: K(z) series is not converged')
        end if
        f(n) = f(n-1)*float(2*n+1)
        t(n) = t(n-1)*u/float(n)
      end do
      do l=1,nmaxlk+6
        f(n+l) = f(n+l-1)*float(2*n+2*l+1)
      end do
      ! Calculation of K(m) with K(0) expansion 
      u = ztab
      do m=1,nmaxlk+6
        ktab(m,i) = u*sum(t(0:n)/f(m:n+m))
        u = u*ztab
      end do
    end do
    ! Recurrence relation factors
    do i=1,nmaxlk+6
      kcoef(i) = float(i)/float(2*i+1)
    end do
  end subroutine

end module
