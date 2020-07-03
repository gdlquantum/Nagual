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
module nmatom
! Manipulation of atomic data
  
  use nmparameter
  use nmtypes
  use nmelement
  use nmshell
  use nmfile
  use nmstring
  use nmunits

  implicit none

    public :: atom_print_geometry
    public :: atom_constructor
    public :: atom_copy
    public :: atom_chk_nelec

    private

    character*(1) :: shell_symbol(10)
    save
    data shell_symbol /'s','p','d','f','g','h','i','j','k','l'/ 

contains

  subroutine atom_print_geometry(tape,atom)
  ! Print geometry for this atom
  ! Roberto Flores-Moreno, 2008
  implicit none
    integer :: tape
    type(natom) :: atom

    character*(2) :: s
    integer :: i

    write (tape,5000) element_get_symbol(atom%atomic_number), &
                      (units_bohr_to_angstrom(atom%pos(i)),i=1,3) 
 5000 format(t2,a2,3x,3f15.6)
  end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Local utilities


  subroutine atom_constructor(this,symbol,pos)
  implicit none
    character*(*) :: symbol
    real(8) :: pos(3)
    type(natom) :: this

    this%symbol(1:2) = symbol(1:2)
    this%pos(1:3) = pos(1:3)
  end subroutine

  subroutine atom_copy(src,dst)
  implicit none
    type(natom) :: src,dst

    dst%symbol(1:2) = src%symbol(1:2)
    dst%pos(1:3) = src%pos(1:3)
    dst%zeff = src%zeff
    dst%atomic_number = src%atomic_number
  end subroutine

  subroutine atom_chk_nelec(this)
  implicit none
    type(natom) :: this

    integer :: nelec

    call element_symbol_to_atomic_number(this%symbol,nelec)
  end subroutine

end module
