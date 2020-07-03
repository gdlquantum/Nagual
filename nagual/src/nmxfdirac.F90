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
module nmxfdirac
! DFT exchange Dirac functional

  use nmmath

  implicit none

    public :: xfdirac_energy
    public :: xfdirac_potential
    public :: xfdirac_kernel
    public :: xfdirac_kernel2

    private 
      real(8) :: r13,r23,r43,r53
      real(8) :: factor0,factor1,factor2,factor3
      parameter ( r13 = 1.0/3.0, r23 = 2.0/3.0, r43 = 4.0/3.0, r53 = 5.0/3.0)
      parameter ( factor0 = -(0.75*(3.0/pi)**r13)*2.0**r43, &
                  factor1 = -((6.0/pi)**r13), &
                  factor2 = -r13*((6.0/pi)**r13), &
                  factor3 = 2.0*r23*((6.0/pi)**r13))

contains
! P. A. M. Dirac, 
! Proc. Cambridge Phil. Soc. 26, 376 (1930).

  subroutine xfdirac_energy(rho,ex)
  ! Evaluates Dirac exchange functional
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    real(8) :: rho,ex

    real(8) :: factor

    ex = factor0*rho**r43
  end subroutine

  subroutine xfdirac_potential(rho,vx)
  ! Evaluates Dirac exchange potential
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    real(8) :: rho,vx

    vx = factor1*rho**r13
  end subroutine

  subroutine xfdirac_kernel(rho,fx)
  ! Evaluates Dirac exchange kernel
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    real(8) :: rho,fx

    integer :: i

    if (rho.eq.0.0) then
      fx = 0.0
    else
      fx = factor2*rho**(-r23)
    end if
  end subroutine

  subroutine xfdirac_kernel2(rho,gx)
  ! Evaluates Dirac exchange second kernel
  ! Roberto Flores-Moreno, 2010, 2018
  implicit none
    real(8) :: rho,gx

    if (rho.eq.0.0) then
      gx = 0.0
    else
      gx = factor3*rho**(-r53)
    end if
  end subroutine

end module 
