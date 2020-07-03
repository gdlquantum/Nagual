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
module nmgui

  use nmtypes
  use nmfile
  use nmsystem
  use nmmolecule
  use nmbasis
  use nmmatrix

  implicit none

    public :: gui_charges
    public :: gui_initialize
    public :: gui_finish
    public :: gui_orbitals
    public :: gui_optimization_step
    public :: gui_perturbation
    public :: gui_perturbation_header

    private
      integer :: ngl

contains

  subroutine gui_initialize(sys)
  ! Initialize interface to GNagual
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys

    call file_open_ascii("Nagual.ngl",ngl)
    call file_copyright(ngl)
    call system_print(ngl,sys)
  end subroutine

  subroutine gui_finish
  ! Finish interface to GNagual
  ! Roberto Flores-Moreno, 2018
  implicit none

    call file_close(ngl)
  end subroutine

  subroutine gui_orbitals(sys)
  ! Print MOs to ngl
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys

    integer :: i,is

    integer :: allocs,dbas,dorb
    real(8),allocatable :: c(:,:)

    write(ngl,'(//,t2,"Orbital energies and coefficients",/)') 
    do is=1,sys%ns
      dbas = basis_nfun(sys%basis(is))
      if (sys%orbital_type.ne.'cartesian') then
        dorb = basis_nfun_spherical(sys%basis(is))
      else
        dorb = dbas
      end if
      write(ngl,'(/,t2,"Set:",i2,/)') is
      write(ngl,'(/,t2,"Number of MOs: ",i5,/)') dorb
      write(ngl,'(/,t2,"Energies for set ",i3,/)') is
      do i=1,dorb
        write(ngl,'(t2,i8,5x,f5.2,5x,f25.13,x,a2)') i,sys%moocc(i,is),&
        sys%moe(i,is)
      end do
      allocate(c(dbas,dbas),stat=allocs)
      if (allocs.gt.0) &
        call file_error('scf/print_orbitals: Allocation failure')
      write(ngl,'(/,t2,"Coefficients for set ",i3,/)') is
      call file_matrix_bas(c,dbas,is,'orbitals','read')
      call matrix_print(ngl,c,dbas,dbas,1,1,dbas,dorb)
      deallocate(c,stat=allocs)
      if (allocs.gt.0) &
        call file_error('scf/print_orbitals: Deallocation failure')
    end do
  end subroutine

  subroutine gui_optimization_step(m,grad,energy,nit)
  ! Print out gradients
  ! Roberto Flores-Moreno, 2013
  implicit none
    integer :: nit
    real(8) :: energy
    real(8) :: grad(3,*)
    type(nmolecule) :: m

    write(ngl,'(/,t2,"Optimization step: ",i5)') nit 
    write(ngl,'(t2,"Energy : ",f20.10)') energy
    call molecule_print_geometry(ngl,m)
  end subroutine 

  subroutine gui_perturbation_header(np)
  ! Print perturbation header
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: np

    write(ngl,'(/,t2,"Responses to density matrices",/)') 
    write(ngl,'(t2,"Number of perturbations: ",i4)') np
  end subroutine
 
  subroutine gui_perturbation(sys,p1,dbas,is,ip)
  ! Print perturbation matrix
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: dbas,is,ip
    real(8) :: p1(dbas,dbas)
    type(nsystem) :: sys

    integer :: i

    write(ngl,'(/,t2,"Density matrix response",/)') 
    write(ngl,'(t2,"Species: ",i4)') is 
    write(ngl,'(t2,"Perturbation: ",i4,/)') ip 
    call matrix_print(ngl,p1,dbas,dbas,1,1,dbas,dbas)
  end subroutine

  subroutine gui_charges(mol,q)
  ! Print charges
  ! Roberto Flores-Moreno, 2019
  implicit none
    real(8) :: q(*)
    type(nmolecule) :: mol

    integer :: iatom

    call file_header(ngl,"Population analysis:") 
    do iatom=1,mol%natom
      write(ngl,'(t2,i3,3x,f14.4)') iatom,q(iatom)
    end do
  end subroutine

end module

