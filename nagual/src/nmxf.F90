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
module nmxf
! Exchange functionals

  use nmparameter
  use nmtypes
  use nmfile
  use nmstring
  use nmxfdirac

  implicit none

    public :: xf_jacob
    public :: xf_code
    public :: xf_energy
    public :: xf_potential
    public :: xf_kernel

    private 

      character*(30) :: xfname(nmaxxf)

      save  ! Fock should be the first 
      data xfname /"FOCK","DIRAC"/

contains

  integer function xf_jacob(sys)
  ! Check if exchange functional is present
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys

    integer :: is,ixf,js
    real(8) :: sex,w

    xf_jacob = 0

    sex = 0.0
    do is=1,sys%ns
      do js=is,sys%ns
        do ixf=2,nmaxxf
          w = abs(sys%interaction(is,js)%exchange(ixf))
          sex = sex + w
          if (w.ne.0.0) then
            xf_jacob = 1
            !if ('pbe96'.eq.string_to_lowercase(xfname(ixf))) xf_jacob = 2
          end if
        end do
      end do
    end do
  end function

  integer function xf_code(str)
  ! Identify exchange functional code (for each new exchange functional
  ! increase nmaxxf (nmparameter) in 1 and add name in xfname)
  ! Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: str

    integer :: ixf

    xf_code = 0
    do ixf=1,nmaxxf
      if (string_to_lowercase(str).eq.&
          string_to_lowercase(xfname(ixf))) xf_code = ixf
    end do
    if (xf_code.eq.0) &
      call file_error('xf_code: unknown exchange functional')
  end function

  subroutine xf_energy(sys,rho,ev,is)
  ! Evaluates exchange functional
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is
    real(8) :: rho,ev
    type(nsystem) :: sys

    integer :: ixf
    real(8) :: t,w

    ev = 0.0

    do ixf=2,nmaxxf   ! 1 for Fock
      w = sys%interaction(is,is)%exchange(ixf)
      if (w.ne.0.0) then
        ! Dirac
        if (xfname(ixf).eq.'DIRAC') then
          call xfdirac_energy(rho,t)
          ev = ev + 0.5*w*t
        !else if (xfname(ixf).eq.'PBE96') then
        !  call xfpbe96_energy(rho,drho,t)
        !  ev = ev + 0.5*w*t
        else if (xfname(ixf).eq.'FOCK') then
          ev = ev + 0.0 ! Added somewhere else
        else
          call file_error('xf_energy: unknown exchange index')
        end if
      end if
    end do
  end subroutine

  subroutine xf_potential(sys,rho,vx,is)
  ! Evaluates exchange potential
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is
    real(8) :: rho,vx
    type(nsystem) :: sys

    integer :: ixf
    real(8) :: t,u,w

    vx = 0.0

    do ixf=2,nmaxxf   ! 1 for Fock
      w = sys%interaction(is,is)%exchange(ixf)
      if (w.ne.0.0) then
        ! Dirac
        if (xfname(ixf).eq.'DIRAC') then
          call xfdirac_potential(rho,t)
          vx = vx + w*t
        !else if (xfname(ixf).eq.'PBE96') then
        !  call xfpbe96_potential(rho,drho,t,u)
        !  ev = ev + 0.5*w*t
        else
          call file_error('xf_potential: unknown exchange index')
        end if
      end if
    end do
  end subroutine

  subroutine xf_kernel(sys,rho,fx,n,is)
  ! Evaluates exchange kernel
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,n
    real(8) :: rho(n),fx(n)
    type(nsystem) :: sys

    integer :: i,ixf
    real(8) :: t,w

    fx(1:n) = 0.0
    
    do ixf=2,nmaxxf   ! 1 for Fock
      w = sys%interaction(is,is)%exchange(ixf)
      if (w.ne.0.0) then
        ! Dirac
        if (ixf.eq.2) then
          do i=1,n
            call xfdirac_kernel(rho(i),t)
            fx(i) = fx(i) + w*t
          end do
        else
          call file_error('xf_kernel: unknown exchange index')
        end if
      end if
    end do
  end subroutine

end module 
