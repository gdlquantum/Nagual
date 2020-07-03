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
module nmopt

  use nmtypes
  use nmfile
  use nmenergy
  use nmlopt
  use nmstring
 
  implicit none

    public :: opt_driver

    private

contains

  subroutine opt_driver(sys,task)
  ! Drive optimization following input directives
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys
    type(ntask) :: task

    integer :: inp_tape

    character*(20) :: algorithm
    integer :: cycles
    namelist /opt/ algorithm,cycles

    call file_open_ascii("Nagual.inp",inp_tape)
    read(inp_tape,nml=opt,end=1000)

    call file_header(task%tape,"Optimization")
    if (string_to_lowercase(algorithm).eq.'local') then
      call file_header(task%tape,"Local relaxation")
      call lopt(cycles,sys,task)
    else
      call file_error('opt_driver, unknown algorithm')
    end if

1000 continue
    call file_close(inp_tape)
  end subroutine

end module 

