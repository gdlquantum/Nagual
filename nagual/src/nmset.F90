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
module nmset
  
  use nmparameter
  use nmtypes
  use nmstring
  use nmfile
  use nmelement
  use nmshell
  use nmunits

  implicit none

    public :: set_lmax
    public :: set_load
    public :: set_get
    public :: set_ecp_get
    public :: set_copy
    public :: set_print

    private

contains

  subroutine set_load(tape,set)
  ! Load one-center basis
  ! Roberto Flores Moreno, 2018
  implicit none
    integer :: tape
    type(nset) :: set  

    character*(1) :: ss
    integer :: i,iblk,n,l,k
    real(8) :: z(nmaxcon),c(nmaxcon)

    read(tape,*) set%nshell
    do iblk=1,set%nshell
      read(tape,*) n,l,k
      if (l.gt.nmaxlb) call file_error('increase nmaxlb')
      if (k.gt.nmaxcon) call file_error('increase nmaxcon')
      do i=1,k
        read(tape,*) z(i),c(i)
      end do
      call shell_load(set%shell(iblk),n,l,k,z,c)
    end do
  end subroutine 

  subroutine set_get(an,set,basfile)
  ! Get atomic basis from library files
  ! Roberto Flores Moreno, 2018
  implicit none
    character*(*) basfile
    integer :: an
    type(nset) :: set

    character*(80) :: line
    character*(32) :: elname
    integer :: tape

    call file_open_ascii(basfile,tape)

    ! Inside file find required basis data
    elname = element_get_name(an)
    read(tape,'(a)',end=99999) line
    do while(string_to_uppercase(line(1:3)).ne.elname(1:3))
      read(tape,'(a)',end=99999) line
    end do

    ! Load basis
    call set_load(tape,set)

    ! Close file
    call file_close(tape)
    return
99999 stop 'basis set not found'
  end subroutine

  subroutine set_ecp_get(an,set,ecpfile,z)
  ! Get atomic ECP from library files
  ! Roberto Flores Moreno, 2020
  implicit none
    character*(*) ecpfile
    integer :: an
    real(8) :: z
    type(nset) :: set

    character*(80) :: line
    character*(32) :: elname
    integer :: ll,ul,tape

    call file_open_ascii(ecpfile,tape)

    ! Inside file find required basis data
    elname = element_get_name(an)
    read(tape,'(a)',end=99999) line
    do while(string_to_uppercase(line(1:3)).ne.elname(1:3))
      read(tape,'(a)',end=99999) line
    end do

    ! Determine number core electrons
    ul = index(line,'|') - 1
    ll = ul - 1
    do while ((line(ll:ll).ge.'0').and.(line(ll:ll).le.'9'))
      ll = ll - 1
    end do
    ll = ll + 1
    read(line(ll:ul),*) z

    ! Load ECPs
    call set_load(tape,set)

    ! Close file
    call file_close(tape)
    return
99999 stop 'ECPs set not found'
  end subroutine

  subroutine set_copy(src,dst)
  ! Copy basis from set src to set dst
  ! Roberto Flores Moreno, 2008, 2018
  implicit none
    type(nset) :: src,dst

    integer :: i

    dst%ll = src%ll
    dst%ul = src%ul
    dst%nshell = src%nshell
    dst%firstrow = src%firstrow
    do i=1,src%nshell
      call shell_copy(src%shell(i),dst%shell(i))
    end do
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Print routines

  subroutine set_print(tape,set)
  ! Print atomic basis to tape
  ! Roberto Flores Moreno, 2018
  implicit none
    integer :: tape
    type(nset) :: set

    integer :: i
    
    write(tape,'(t2,"===============")') 
    write(tape,*) set%nshell
    do i=1,set%nshell
      call shell_print(tape,set%shell(i))
    end do
  end subroutine 

  integer function set_lmax(set)
  ! Obtain maximum l for a given set
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nset) :: set

    integer :: ishell,l

    set_lmax = 0
    do ishell=1,set%nshell
      l = set%shell(ishell)%l
      if (l.gt.set_lmax) set_lmax = l
    end do
  end function

end module
