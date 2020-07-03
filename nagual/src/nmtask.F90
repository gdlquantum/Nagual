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
module nmtask
! Task information
  
  use nmparameter
  use nmtypes
  use nmfile
  use nmstring
  use nmtime
  use nmmath
  use nmintegrals
  use nmscf

  use nmsystem
  use nmbasis
  use nmmolecule

  use nmenergy
  use nmgreen
  use nmgui
  use nmaux
  use nmopt

  implicit none

    public :: task_run
    public :: task_load_from_file

    private

contains

  subroutine task_run(this,sys)
  ! Run task 
  ! Roberto Flores-Moreno, 2008, 2009, 2014, 2018
  implicit none
    type(ntask) :: this
    type(nsystem) :: sys

    integer :: is

    ! Initialize 
    call integrals_initialize(sys%mol,this%tape)
    do is=1,sys%ns
      call basis_initialize(sys%basis(is))
    end do
    if (this%scf%eris.eq.3) call aux_initialize(this%tape,sys)

    ! Initial infor to GUI file
    call gui_initialize(sys)

    ! Geometry optimization
    call opt_driver(sys,this)

    ! Energy
    call energy_single_point(this%tape,this,sys,sys%energy,.true.)
    call gui_orbitals(sys)

    call green_driver(this%tape,sys)

    call gui_finish
  end subroutine

  subroutine task_load_from_file(tape,task,sys)
  ! Get task data from input ascii file
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape
    type(ntask) :: task
    type(nsystem) :: sys

    character*(80) :: line
    integer :: i,n

    character*(20) :: guess
    integer :: diis,eris,nitmax
    real(8) :: damping,tolerance
    namelist /scf/ diis,nitmax,guess,eris,damping,tolerance

    integer :: species_id,orbital_id
    namelist /tom/ species_id,orbital_id

    character*(20) :: family
    integer :: order
    namelist /sq/ family,order

    task%sq%family = 'skip'
    task%sq%order = 0
    rewind(tape)
    read(tape,NML=sq,end=1000)
    task%sq%family = family
    task%sq%order = order
    task%scf%guess = 'skip'
    return
 1000 continue

    ! Load molecular system if required
    call system_load_from_file(tape,task%tape,sys)

    ! TOM (Slater's transition operator method)
    species_id = 0
    orbital_id = 0
    rewind(tape)
    read(tape,NML=tom,end=1100)
 1100 continue
    if ((species_id.gt.0).and.(orbital_id.gt.0)) then
      if (species_id.le.sys%ns) then
        if (orbital_id.le.sys%nmoo(species_id)) then
          sys%moocc(orbital_id,species_id) = 0.5
        end if
      end if
    end if

    ! Guess MOs
    task%scf%nitmax = 100
    task%scf%guess = 'core'  ! Default to guess core
    task%scf%eris = 4 ! Default to exact 4-ERIs
    task%scf%damping = 1.0
    task%scf%tol = 1.0e-8
    task%scf%diis = -1
    rewind(tape)
    read(tape,nml=scf,end=1300) 
    task%scf%guess = string_to_lowercase(guess)
    task%scf%eris = eris
    task%scf%nitmax = nitmax
    task%scf%tol = tolerance
    task%scf%diis = diis
    if (diis.gt.0) then
      task%scf%damping = 1.0
    else
      task%scf%damping = damping
    end if
 1300 continue
  end subroutine

end module
