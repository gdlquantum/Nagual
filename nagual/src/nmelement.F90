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
module nmelement
  
  use nmparameter
  use nmfile
  use nmstring
  use nmunits

  implicit none
   
    public :: element_get_symbol
    public :: element_get_name
    public :: element_symbol_to_atomic_number
    public :: element_mass_from_atomic_number

    private

      character*(2) :: element_symbol(0:nmaxatnum)
      character*(30) :: element_name(0:nmaxatnum)
      real(8) :: element_mass(0:nmaxatnum)

      save
      data element_symbol /'X',&
      'H','He',&
      'Li','Be',&
      'B','C','N','O','F','Ne',&
      'Na','Mg',&
      'Al','Si','P','S','Cl','Ar',&
      'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
      'Ga','Ge','As','Se','Br','Kr',&
      'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',&
      'In','Sn','Sb','Te','I','Xe',&
      'Cs','Ba',& 
      'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',&
      'Hf','Ta',' W','Re','Os','Ir','Pt','Au',&
      'Hg','Tl','Pb','Bi','Po','At','Rn',&
      'Fr','Ra',&
      'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/

      ! Lit.: CRC Handbook of Chemistry and Physics, 1989 
      data element_mass  &
                /  0.000000 ,  1.007940 ,  4.002602 ,  6.941000 ,  9.012182&
                , 10.811000 , 12.011000 , 14.006740 , 15.999400 , 18.998400&
                , 20.179700 , 22.989768 , 24.305000 , 26.981539 , 28.085500&
                , 30.973762 , 32.066000 , 35.452700 , 39.948000 , 39.098300&
                , 40.078000 , 44.955910 , 47.880000 , 50.941500 , 51.996100&
                , 54.938050 , 55.847000 , 58.933200 , 58.693400 , 63.546000&
                , 65.390000 , 69.723000 , 72.610000 , 74.921590 , 78.960000&
                , 79.904000 , 83.800000 , 85.467800 , 87.620000 , 88.905850&
                , 91.224000 , 92.906380 , 95.940000 , 98.000000 ,101.070000&
                ,102.905500 ,106.420000 ,107.868200 ,112.411000 ,114.820000&
                ,118.710000 ,121.757000 ,127.600000 ,126.904470 ,131.290000&
                ,132.905430 ,137.327000 ,138.905500 ,140.115000 ,140.907650&
                ,144.240000 ,145.000000 ,150.360000 ,151.965000 ,157.250000&
                ,158.925340 ,162.500000 ,164.930320 ,167.260000 ,168.934210&
                ,173.040000 ,174.967000 ,178.490000 ,180.947900 ,183.850000&
                ,186.207000 ,190.200000 ,192.220000 ,195.080000 ,196.966540&
                ,200.590000 ,204.383300 ,207.200000 ,208.980370 ,209.000000&
                ,210.000000 ,222.000000 ,223.000000 ,226.000000 ,227.000000&
                ,232.038100 ,231.035880 ,238.028900 ,237.000000 ,244.000000&
                ,243.000000 ,247.000000 ,247.000000 ,251.000000 ,252.000000&
                ,257.000000 ,258.000000 ,259.000000 ,262.000000/

      data element_name  &
                /"DUMMY"      ,"HYDROGEN"      ,"HELLIUM"    ,"LITHIUM"     &
                ,"BERYLLIUM"  ,"BORON"         ,"CARBON"     ,"NITROGEN"    &
                ,"OXYGEN"     ,"FLUORINE"      ,"NEON"       ,"SODIUM"      &
                ,"MAGNESIUM"  ,"ALUMINIUM"     ,"SILICON"    ,"PHOSPHORUS"  &
                ,"SULFUR"     ,"CHLORINE"      ,"ARGON"      ,"POTASSIUM"   &
                ,"CALCIUM"    ,"SCANDIUM"      ,"TITANIUM"   ,"VANADIUM"    &
                ,"CHROMIUM"   ,"MANGANESE"     ,"IRON"       ,"COBALT"      &
                ,"NICKEL"     ,"COPPER"        ,"ZINC"       ,"GALLIUM"     &
                ,"GERMANIUM"  ,"ARSENIC"       ,"SELENIUM"   ,"BROMINE"     &
                ,"KRYPTON"    ,"RUBIDIUM"      ,"STRONTIUM"  ,"YTTRIUM"     &
                ,"ZIRCONIUM"  ,"NIOBIUM"       ,"MOLYBDENUM" ,"TECHNETIUM"  &
                ,"RUTHENIUM"  ,"RHODIUM"       ,"PALLADIUM"  ,"SILVER"      &
                ,"CADMIUM"    ,"INDIUM"        ,"TIN"        ,"ANTIMONY"    &
                ,"TELLURIUM"  ,"IODINE"        ,"XENON"      ,"CESIUM"      &
                ,"BARIUM"     ,"LANTHANUM"     ,"CERIUM"     ,"PRASEODYMIUM"&
                ,"NEODYMIUM"  ,"PROMETHIUM"    ,"SAMARIUM"   ,"EUROPIUM"    &
                ,"GADOLINIUM" ,"TERBIUM"       ,"DYSPROSIUM" ,"HOLMIUM"     &
                ,"ERBIUM"     ,"THULIUM"       ,"YTTERBIUM"  ,"LUTETIUM"    &
                ,"HAFNIUM"    ,"TANTALUM"      ,"TUNGSTEN"   ,"RHENIUM"     &
                ,"OSMIUM"     ,"IRIDIUM"       ,"PLATINUM"   ,"GOLD"        &
                ,"MERCURY"    ,"THALLIUM"      ,"LEAD"       ,"BISMUTH"     &
                ,"POLONIUM"   ,"ASTATINE"      ,"RADON"      ,"FRANCIUM"    &
                ,"RADIUM"     ,"ACTINIUM"      ,"THORIUM"    ,"PROTACTINIUM"&
                ,"URANIUM"    ,"NEPTUNIUM"     ,"PLUTONIUM"  ,"AMERICIUM"   &
                ,"CURIUM"     ,"BERKELIUM"     ,"CALIFORNIUM","EINSTEINIUM" &
                ,"FERMIUM"    ,"MENDELEVIUM"   ,"NOBELIUM"   ,"LAWRENCIUM"  /

contains

  character*(2) function element_get_symbol(atomic_number)
  ! Access to symbol
  ! Roberto Flores-Moreno (Oct 2014)
  implicit none
    integer :: atomic_number

    element_get_symbol(1:2) = element_symbol(atomic_number)(1:2)
  end function

  character*(30) function element_get_name(atomic_number)
  ! Access to symbol
  ! Roberto Flores-Moreno 2018
  implicit none
    integer :: atomic_number

    element_get_name(1:30) = element_name(atomic_number)(1:30)
  end function

  subroutine element_symbol_to_atomic_number(symbol,atomic_number)
  ! Access to symbol if atomic number is available
  ! Roberto Flores-Moreno (Oct 2014)
  implicit none
    character*(2) :: symbol
    integer :: atomic_number

    integer :: i

    atomic_number = 0
    do i=0,nmaxatnum
      if (string_to_lowercase(symbol).eq.&
          string_to_lowercase(element_symbol(i))) then
        atomic_number = i
      end if
    end do
  end subroutine

  subroutine element_mass_from_atomic_number(atomic_number,mass)
  ! Access to mass
  ! Roberto Flores-Moreno,  2019
  implicit none
    integer :: atomic_number
    real(8) :: mass

    integer :: i

    mass = 0.0

    do i=0,nmaxatnum
      if (atomic_number.eq.i) then 
        mass = element_mass(i)
      end if
    end do
  end subroutine

end module
