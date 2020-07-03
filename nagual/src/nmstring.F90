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
module nmstring
! Management of character strings

  implicit none

    public :: string_to_lowercase
    public :: string_to_uppercase
    public :: string_size
    public :: string_from_integer

    private

contains

  character*(256) function string_to_lowercase(string)
  ! Purpose: Converts all upper case characters in a string to
  !          lower case characters.
  ! Roberto Flores Moreno, Jul 2008, Sep 2008
  implicit none
    character :: string*(*)
    integer :: i,iascii
    
    string_to_lowercase = string
    do i=1,len(string)
      iascii = ichar(string(i:i))
      if ((iascii>=65).and.(iascii<=90)) string_to_lowercase(i:i) = char(iascii+32)
    end do
  end function string_to_lowercase

  character*(256) function string_to_uppercase(string)
  ! Purpose: Converts all lower case characters in a string to
  !          upper case characters.
  ! Roberto Flores Moreno, Sep 2008
  implicit none
    character :: string*(*)
    integer :: i,iascii

    string_to_uppercase = string
    do i=1,len(string)
      iascii = ichar(string(i:i))
      if ((iascii.ge.97).and.(iascii.le.122)) &
        string_to_uppercase(i:i) = char(iascii - 32)
    end do
  end function

  integer function string_size(string)
  ! Purpose: Get STRing EXTension.
  ! Roberto Flores Moreno 2008
  implicit none
    character :: string*(*)
    integer :: i

    string_size = 1
    do i=len(string),1,-1
      if (string(i:i).ne.' ') then
        string_size = i
        return
      end if
    end do
  end function 

  subroutine string_from_integer(n,str)
  ! Build compact string for an integer number
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: n
    character*(*) :: str
    
    integer :: m
    character*(50) :: fs

    if (n.lt.0) STOP 'Negative number in string_from_integer'

    m = int(log10(float(n))) + 1
    write(fs,'(t1,"(t1,i",i1,")")') m
    write(str,fs) n
  end subroutine

end module
