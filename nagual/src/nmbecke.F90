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
module nmbecke
! Handling of molecular grids

  use nmparameter
  use nmtypes
  use nmmath
  use nmvector

  implicit none

    public :: becke_radial_quadrature
    public :: becke_skgc
    public :: becke_tskgc
    public :: becke_weight

    private 

      integer :: nradp,nangp
      parameter (nradp = 75,& ! Number of radial points / atom
                 nangp = 590)  ! Number of angular points / shell

      real(8) :: rg(4,nradp*nangp)
      real(8) :: g(4,nradp*nangp)

contains

  real(8) function becke_function(mu)
  ! Becke iteration function
  ! Roberto Flores-Moreno, 2010
  implicit none
    real(8) :: mu

    becke_function = 0.5*mu*(3.0 - mu**2)
  end function

  subroutine becke_radial_quadrature(r,w,n)
  ! Radial quadrature
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    integer :: n
    real(8) :: r(*),w(*)

    ! Transformed Gauss-Chebyshev quadratures
    call becke_tskgc(r,w,n)

    w(:n) = (w(:n)/(1.0-r(:n)))/log(2.0)
    r(:n) = 1.0 - log(1.0-r(:n))/log(2.0)
  end subroutine

  subroutine becke_skgc(x,w,n)
  ! Transformed second kind Gauss-Chebyshev quadratures
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: n
    real(8) :: x(*),w(*)
  
    integer :: i
    real(8) :: ct,st,t

    do i=1,n
      t = pi*float(i)/float(n+1)
      ct = cos(t)
      st = sin(t)
      x(i) = ct
      w(i) = st*t/float(i) 
    end do
  end subroutine

  subroutine becke_tskgc(x,w,n)
  ! Transformed second kind Gauss-Chebyshev quadratures
  ! Roberto Flores-Moreno, May 2009
  !
  ! J. M. Perez-Jorda, E. San-Fabian, F. Moscardo,
  ! Comput. Phys. Commun. 70, 271 (1992).
  !
  ! J. M. Perez-Jorda, A. Becke, E. San-Fabian,
  ! J. Chem. Phys. 100, 6520 (1994).
  implicit none
    integer :: n
    real(8) :: x(*),w(*)
  
    integer :: i
    real(8) :: ct,st,st2,t

    do i=1,int((n+1)/2)
      t = pi*float(i)/float(n+1)
      ct = cos(t)
      st = sin(t)
      st2 = st*st
      x(i) = float(2*i-n-1)/float(n+1)-2.0*(1.0+2.0*st2/3.0)*ct*st/pi
      w(i) = 16.0*st2**2/float(3*(n+1))
      x(n-i+1) = -x(i)
      w(n-i+1) = w(i)
    end do
  end subroutine

  real(8) function becke_weight(m,atoma,pos)
  ! Evaluates Becke weight for a given point
  ! A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    integer :: atoma
    real(8) :: pos(3)

    integer :: iatom,jatom
    real(8) :: mu
    real(8) :: r1(3),r2(3)
    real(8) :: r(nmaxatom),p(nmaxatom)

    do iatom=1,m%natom
      r1(1:3) = m%atom(iatom)%pos(1:3)
      r(iatom) = vector_distance(pos,r1)
      p(iatom) = 1.0
    end do
    do iatom=2,m%natom
      r1(1:3) = m%atom(iatom)%pos(1:3)
      do jatom=1,iatom-1
        r2(1:3) = m%atom(jatom)%pos(1:3)
        mu = vector_distance(r2,r1)
        mu = (r(iatom)-r(jatom))/mu
        mu = becke_function(becke_function(becke_function(mu)))
        p(iatom) = (1.0 - mu)*p(iatom)
        p(jatom) = (1.0 + mu)*p(jatom)
      end do
    end do
    becke_weight = p(atoma)/sum(p(1:m%natom))
  end function

end module 
