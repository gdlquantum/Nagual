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
module nmmath
! Some mathematical functions and utilities
!

  implicit none

    public :: math_factorial
    public :: math_factoriald
    public :: math_binomial

    real(8), public :: pi
    parameter ( pi = 2.0*acos(0.0) )

    private

contains

  real(8) function math_binomial(n,k)
  ! Evaluate binomial coefficient
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: k,n

    math_binomial = math_factorial(n)/(math_factorial(n-k)*math_factorial(k))
  end function

  real(8) function math_factorial(n)
  ! Evaluate factorial n!
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: n
  
    integer :: i

    if (n.lt.0) then
      math_factorial = 0.0
      return
    end if
    math_factorial = 1.0
    do i=2,n
      math_factorial = math_factorial*float(i)
    end do
  end function 

  real(8) function math_factoriald(n)
  ! Evaluate double factorial n!!
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: n
  
    integer :: i

    if (n<-1) then
      math_factoriald = 0.0
      return
    end if
    math_factoriald = 1.0
    if (n==2) math_factoriald = 2.0
    do i=mod(n,2)+2,n,2
      math_factoriald = math_factoriald*float(i)
    end do
  end function 

end module
