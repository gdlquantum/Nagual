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
module nmtypes

  use nmparameter

  implicit none

    type, public :: ntimer
      character*(50) :: title
      integer :: stat                    ! 0: off 1: on
      real(4) :: val
      real(4) :: when
    end type ntimer

    type, public :: nshell                    ! A basis set shell
      integer :: n                            ! level/radial power
      integer :: l                            ! Angularity
      integer :: k                            ! Contraction degree
      real(8) :: radius
      real(8) :: z(nmaxcon)                   ! Gaussian exponents
      real(8) :: c(nmaxcon)                   ! Contraction coefficients
      real(8) :: norm(nmaxncok)                ! Normalization for AOs
    end type nshell

    ! Each set is located on atom with same number
    type, public :: nset
      integer :: atom
      integer :: ll 
      integer :: ul
      integer :: nshell                       ! Number of shells
      logical :: firstrow
      real(8) :: radius
      type(nshell) :: shell(nmaxashell)        ! Atoms
    end type nset

    type, public :: nbasis
      character*(256) :: filename
      integer :: nsets
      type(nset) :: set(nmaxset)      
    end type nbasis

    type, public :: natom ! Atom structure to collect data
      character*(2) :: symbol
      integer :: atomic_number             ! Type
      real(8) :: zeff                    ! Effective zet 
      integer :: zref(3)                    
      real(8) :: zpos(3)      
      real(8) :: pos(3)                    ! Location
    end type natom

    type, public :: nmolecule   ! Molecule object 
      character*(3) :: pg
      integer :: natom                        ! Number of atoms
      type(natom) :: atom(nmaxatom)                     ! Atoms
    end type nmolecule

    ! Order for interaction types 
    ! simple: (cc,yuk)
    !          1) cc: Coulomb coefficient 
    !          2) yuk : Yukawa amplifier
    type, public :: ninteraction
      real(8) :: simple(2) 
      real(8) :: exchange(nmaxxf)
      real(8) :: correlation(nmaxcf)
    end type ninteraction

    type, public :: nsystem   
      character*(20) :: orbital_type             ! Spherical or cartesian
      character*(20) :: species_name(nmaxns)  ! Particle names
      type(ninteraction) :: interaction(nmaxns,nmaxns)
      type(nmolecule) :: mol                    ! Molecule structure and basis
      character*(20) :: species(nmaxns)  
      integer :: ns                             ! Number of species types
      integer :: nmoo(nmaxns)                  ! Occupied MOs
      real(8) :: species_charge(nmaxns)       ! Particle charges
      real(8) :: species_mass(nmaxns)         ! Particle masses
      real(8) :: energy                          
      real(8) :: moe(nmaxnbas,nmaxns)          ! MO energies
      real(8) :: moocc(nmaxnbas,nmaxns)        ! MO occupations
      type(nbasis) :: basis(nmaxns)
      type(nbasis) :: auxis(nmaxns)
      type(nbasis) :: ecps
    end type nsystem

    type, public :: nscftask                    ! SCF task description
      character*(20) :: guess                   ! Guess type for SCF
      integer :: diis                           ! Number of matrices in DIIS
      integer :: eris                           ! 4-ERIs or 3-ERIs
      integer :: nitmax
      real(8) :: damping                        ! damping value for scf
      real(8) :: tol
    end type nscftask

    type, public :: nsqtask                    ! SCF task description
      character*(20) :: family
      integer :: order
    end type nsqtask

    type, public :: ntask                       ! Job description
      integer :: tape                              ! Print tape
      type(nscftask) :: scf                        ! SCF 
      type(nsqtask) :: sq                          ! Second Quantization
    end type ntask

end module
