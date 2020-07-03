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
module nmgreen

  use nmtypes
  use nmfile
  use nmg1pp

  implicit none

    public :: green_driver

    private

contains

  subroutine green_driver(tape,sys)
  ! Driver for Green's function calculations
  ! Roberto Flores-Moreno 2015, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys
   
    call g1pp_wraper(tape,sys)
  end subroutine 

  subroutine g1pp_wraper(tape,sys)
  ! Driver for G1PP calculations
  ! Roberto Flores-Moreno 2014, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: ist,inp_tape

    character*(20) :: method
    integer :: memory,species_id,orbital_id
    namelist /g1pp/ species_id,orbital_id,method,memory

    species_id = 0
    call file_open_ascii("Nagual.inp",inp_tape)
    read(inp_tape,nml=g1pp,end=1000)
1000 continue
    call file_close(inp_tape)
    if (species_id.le.0) return

    call g1pp_driver(tape,sys,species_id,orbital_id,method,memory)
  end subroutine

end module
