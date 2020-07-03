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
module nmparameter
! Define parameters 

  implicit none

    public 

    integer :: nmaxatnum             ! maximum atomic number allowed
    integer :: nmaxatom              ! maximum number of atoms
    integer :: nmaxcon               ! maximum degree of contraction for basis
    integer :: nmaxcf                ! maximum number of correlation functionals
    integer :: nmaxgpa               ! maximum grid points per atom
    integer :: nmaxl                 ! maximum angularity for integration
    integer :: nmaxlb                ! maximum l for basis
    integer :: nmaxlp                ! maximum l for potential
    integer :: nmaxli                ! maximum angularity (composed)
    integer :: nmaxlk                ! maximum angularity for bessel functions
    integer :: nmaxlder              ! maximum order of derivatives
    integer :: nmaxnbas              ! maximum number of basis functions
    integer :: nmaxnchem             ! maximum number of chemicals
    integer :: nmaxnco               ! maximum number of basis in a shell
    integer :: nmaxncok        
    integer :: nmaxnfio              ! maximum number of files for I/O
    integer :: nmaxnisotopes         ! maximum number of isotopes per element
    integer :: nmaxnreact            ! maximum number of reactants in reaction
    integer :: nmaxnrxns             ! maximum number of reactions
    integer :: nmaxns                ! maximum number of species types
    integer :: nmaxram               ! maximum RAM in megabytes
    integer :: nmaxset               ! maximum number of sets per basis
    integer :: nmaxshell             ! maximum number of basis shells 
    integer :: nmaxashell            ! maximum number of basis shells in set
    integer :: nmaxtgam              ! maximum number of tabulated gamma values
    integer :: nmaxxf                ! maximum number of exchange functionals

    real(8) :: ntolnum               ! tiniest non negligible number
    real(8) :: ntolint               ! Integration tolerance
    real(8) :: ntolgeo               ! Geometry tolerance

    ! Physical constants
    real(8) :: clight
    real(8) :: hplanck
    real(8) :: kboltzmann
    real(8) :: navogadro

    parameter (nmaxatnum = 103)
    parameter (nmaxnisotopes = 6)
    parameter (nmaxnfio = 20)

    parameter (nmaxatom = 100)
    parameter (nmaxcon = 25)
    parameter (nmaxcf = 1)
    parameter (nmaxgpa = 30000)
    parameter (nmaxlb = 4)             
    parameter (nmaxlp = 5)             
    parameter (nmaxlder = 0)
    parameter (nmaxtgam = 120)
    parameter (nmaxset = nmaxatom)
    parameter (nmaxashell = 50)
    parameter (nmaxns = 5)
    parameter (nmaxl = nmaxlb + nmaxlder)
    parameter (nmaxli = nmaxl + max(nmaxlb,nmaxlp))  
    parameter (nmaxlk = nmaxli + nmaxlp)
    parameter (nmaxnco = ((nmaxli+1)*(nmaxli+2))/2)
    parameter (nmaxncok = ((nmaxli+3)*(nmaxli+4))/2)
    parameter (nmaxshell = nmaxset*nmaxashell)
    parameter (nmaxnbas = nmaxnco*nmaxshell)
    parameter (nmaxxf = 2)

    parameter (ntolnum = 1.0e-14)
    parameter (ntolint = ntolnum)
    parameter (ntolgeo = 1.0e-4)

    parameter (navogadro = 6.0221367e+23, & 
               clight = 299792458.0, &       
               hplanck = 6.62606957e-34*navogadro, &    ! J/mol/s
               kboltzmann = 8.3144621)              ! J/mol/K
end module 
