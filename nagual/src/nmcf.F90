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
module nmcf
! Correlation functionals

  use nmparameter
  use nmtypes
  use nmfile
  use nmstring
  use nmcfvwn

  implicit none

    public :: cf_jacob
    public :: cf_code
    public :: cf_energy
    public :: cf_potential
    public :: cf_kernel

    private 

      character*(30) :: cfname(nmaxcf)

      save
      data cfname /"VWN"/

contains

  integer function cf_jacob(sys)
  ! Check if correlation functional is present
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys

    integer :: is,icf,js
    real(8) :: sec,w

    cf_jacob = 0

    sec = 0.0
    do is=1,sys%ns
      do js=is,sys%ns
        do icf=1,nmaxcf
          w = abs(sys%interaction(is,js)%correlation(icf))
          sec = sec + w
          if (w.ne.0.0) then
            cf_jacob = 1
          end if
        end do
      end do
    end do
  end function

  integer function cf_code(str)
  ! Identify correlation functional code (for each new correlation functional
  ! increase nmaxcf (nmparameter) in 1 and add name in cfname)
  ! Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: str

    integer :: icf

    cf_code = 0
    do icf=1,nmaxcf
      if (string_to_lowercase(str).eq.&
          string_to_lowercase(cfname(icf))) cf_code = icf
    end do
    if (cf_code.eq.0) &
      call file_error('cf_code: unknown correlation functional')
  end function

  subroutine cf_energy(sys,rhoa,rhob,ev,is,js)
  ! Evaluates correlation energy functional
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,js
    real(8) :: rhoa,rhob,ev
    type(nsystem) :: sys

    integer :: icf
    real(8) :: t,w

    ev = 0.0

    do icf=1,nmaxcf
      w = sys%interaction(is,js)%correlation(icf)
      if (w.ne.0.0) then
        if (cfname(icf).eq.'VWN') then 
          call cfvwn_energy(rhoa,rhob,t)
          ev = ev + w*t
        else
          call file_error('cf_energy: unknown functional index')
        end if
      end if
    end do
  end subroutine

  subroutine cf_potential(sys,rhoa,rhob,vca,vcb,is,js)
  ! Evaluates correlation potential
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,js
    real(8) :: rhoa,rhob,vca,vcb
    type(nsystem) :: sys

    integer :: icf
    real(8) :: ta,tb,w

    vca = 0.0
    vcb = 0.0

    do icf=1,nmaxcf
      w = sys%interaction(is,js)%correlation(icf)
      if (w.ne.0.0) then
        if (icf.eq.1) then 
          call cfvwn_potential(rhoa,rhob,ta,tb)
          vca = vca + w*ta
          vcb = vcb + w*tb
        else
          call file_error('cf_potential: unknown functional index')
        end if
      end if
    end do
  end subroutine

  subroutine cf_kernel(sys,rhoa,rhob,fcaa,fcab,fcbb,n,is,js)
  ! Evaluates correlation kernel
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,js,n
    real(8) :: rhoa(n),rhob(n),fcaa(n),fcab(n),fcbb(n)
    type(nsystem) :: sys

    integer :: i,icf
    real(8) :: taa,tab,tbb,w

    fcaa(1:n) = 0.0
    fcab(1:n) = 0.0
    fcbb(1:n) = 0.0

    do icf=1,nmaxcf
      w = sys%interaction(is,js)%correlation(icf)
      if (w.ne.0.0) then
        !if (icf.eq.1) then
        !  do i=1,n
        !    call cfvwn_kernel(rhoa(i),rhob(i),taa,tab,tbb)
        !    fcaa(i) = fcaa(i) + w*taa
        !    fcab(i) = fcab(i) + w*tab
        !    fcbb(i) = fcbb(i) + w*tbb
        !  end do
        !else
          call file_error('cf_kernel: unknown functional index')
        !end if
      end if
    end do
  end subroutine

end module 
