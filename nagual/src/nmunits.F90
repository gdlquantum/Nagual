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
module nmunits
! Units conversion factors

  use nmparameter

  implicit none

    public :: units_jmol_to_au
    public :: units_au_to_jmol
    public :: units_au_to_kcalmol
    public :: units_wavenumber_to_jmol
    public :: units_angstrom_to_bohr
    public :: units_bohr_to_angstrom
    public :: units_au_to_ev
    public :: units_au_to_esu
    public :: units_au_to_debye
    public :: units_ev_to_au

    private

contains

  real(8) function units_au_to_ev(au)
  ! Transform au to electron volts
  ! Roberto Flores Moreno, 2008
  implicit none
    real(8) :: au

    units_au_to_ev = 27.21138*au
  end function 

  real(8) function units_ev_to_au(ev)
  ! Transform electron volts to au
  ! Roberto Flores Moreno, 2018
  implicit none
    real(8) :: ev

    units_ev_to_au = ev/27.21138
  end function 

  real(8) function units_angstrom_to_bohr(aa)
  ! Transform angstrom to bohr
  ! Roberto Flores Moreno, 2008
  implicit none
    real(8) :: aa

    units_angstrom_to_bohr = aa/0.5291772080172567
  end function 

  real(8) function units_bohr_to_angstrom(db)
  ! Transform angstrom to bohr
  ! Roberto Flores Moreno, (Oct 2014)
  implicit none
    real(8) :: db

    units_bohr_to_angstrom = db*0.5291772080172567
  end function 

  real(8) function units_au_to_jmol(au)
  ! Transform au to J/mol
  ! Roberto Flores-Moreno 2018
  implicit none
    real(8) :: au

    units_au_to_jmol = au*navogadro*4.35974417e-18
  end function 

  real(8) function units_jmol_to_au(jmol)
  ! Transform J/mol to au
  ! Roberto Flores-Moreno 2018
  implicit none
    real(8) :: jmol

    units_jmol_to_au = jmol/(navogadro*4.35974417e-18)
  end function 

  real(8) function units_au_to_kcalmol(au)
  ! Transform au to J/mol
  ! Roberto Flores-Moreno 2018
  implicit none
    real(8) :: au

    units_au_to_kcalmol = au*627.5 !FIXME
  end function 

  real(8) function units_wavenumber_to_jmol(wn)
  ! Transform cm^-1 to J/mol
  implicit none
    real(8) :: wn

    units_wavenumber_to_jmol = 100.0*wn*hplanck*clight
  end function 

  real(8) function units_au_to_esu(au)
  ! Transform au to ESU
  ! Roberto Flores Moreno, 2019
  implicit none
    real(8) :: au

    units_au_to_esu = units_bohr_to_angstrom(au*3.0*1.6021765314)
  end function 

  real(8) function units_au_to_debye(au)
  ! Transform au to Debye
  ! Roberto Flores Moreno, 2020
  implicit none
    real(8) :: au

    units_au_to_debye = au/0.3935
  end function 

end module
