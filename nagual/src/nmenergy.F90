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
module nmenergy

  use nmtypes
  use nmfile
  use nmunits
  use nmtime
  use nmscf
  use nmgrid
  use nmxf
  use nmcf
  use nmmp

  implicit none

    public :: energy_single_point

    private

contains

  subroutine energy_single_point(tape,task,sys,e,report)
  ! Compute energy for a given structure with different methods
  ! Roberto Flores-Moreno, 2014, 2018
  implicit none
    integer :: tape
    logical :: report
    real(8) :: e
    type(nsystem) :: sys
    type(ntask) :: task

    logical :: dft
    type(ntimer) :: tused

    if (report) then
      call start_timer(tused,'Energy time:')
      write(tape,'(t2,"=== Energy calculation ===")') 
    end if

    dft = .false.
    if (xf_jacob(sys)+cf_jacob(sys).gt.0) dft = .true.
    if (dft) call grid_initialize(sys)

    if (.not.scf(task%scf,sys)) then
      call file_error('energy_sigle_point, SCF did not converged')
    else if (report) then
      write(tape,'(t2,"SCF ENERGY (a.u.): ",f30.6)') sys%energy
      write(tape,'(t2,"             (eV): ",f29.5)') &
        units_au_to_ev(sys%energy)
      write(tape,'(t2,"       (Kcal/mol): ",f27.3)') &
        units_au_to_kcalmol(sys%energy)
      write(tape,'(t2,"For more details see SCF file",/)') 
    end if

    ! Moller-Plesset MBPT
    call mp_driver(tape,sys)
    if (report) write(tape,'(t2,"TOTAL ENERGY (a.u.): ",f30.6)') sys%energy

    e = sys%energy

    if (report) call print_timer(tused,tape)
  end subroutine

end module
