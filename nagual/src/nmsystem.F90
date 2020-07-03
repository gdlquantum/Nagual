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
module nmsystem
! System administration and definitions
  
  use nmparameter
  use nmtypes
  use nmfile
  use nmstring
  use nmunits
  use nmelement
  use nmatom
  use nmmolecule
  use nmbasis
  use nmset
  use nmxf
  use nmcf

  implicit none

    public :: system_load_from_file
    public :: system_multiplicity
    public :: system_charge
    public :: system_print

    private

contains

  subroutine system_load_from_file(tape,outf,sys)
  ! Get system data from input ascii file
  ! Roberto Flores Moreno, Jul 2008, Nov 2008
  implicit none
    integer :: tape,outf
    type(nsystem) :: sys

    character*(80) :: line,str,pg
    character*(20) :: name1,name2,name3
    character*(10) :: symbol
    integer :: i,iatom,is,n,n2,nbas,nelec,l,unpaired
    real(8) :: val
    logical :: inbohr
    type(nmolecule) :: m
    type(natom) :: atom
    type(nset) :: set

    character*(254) :: basfile,ecpfile

    character*(20) :: classical,exchange,correlation

    !!! Order is very important in reading 

    ! Geometry 
    rewind(tape)
    read(tape,'(a)',end=99999) line
    do while(string_to_lowercase(line(1:10)).ne.'&structure')
      read(tape,'(a)',end=99999) line
    end do
    if (index(string_to_lowercase(line),'bohr').ne.0) then
      inbohr = .true.
    else
      inbohr = .false.
    end if
    if (index(string_to_lowercase(line),'z').ne.0) then
      call molecule_read_zmatrix(sys%mol,tape,'/',inbohr)
    else
      call molecule_read_cartesian(sys%mol,tape,'/',inbohr)
    end if

    ! ECPs
    rewind(tape)
    read(tape,'(a)',end=1200) line
    do while(string_to_lowercase(line(1:5)).ne.'&ecps')
      read(tape,'(a)',end=1200) line
    end do
    read(tape,*) ecpfile
    call get_ecps(sys,ecpfile)
 1200 continue

    ! Species
    rewind(tape)
    read(tape,'(a)',end=1300) line
    do while(string_to_lowercase(line(1:10)).ne.'&species')
      read(tape,'(a)',end=1300) line
    end do
    read(tape,*) sys%ns
    do i=1,sys%ns
      read(tape,*) sys%species_name(i),sys%nmoo(i),basfile
      if (string_to_lowercase(sys%species_name(i)).eq.'electron') then
        sys%species_charge(i) = -1.0
        sys%species_mass(i) = 1.0
      else if (string_to_lowercase(sys%species_name(i)).eq.'proton') then
        sys%species_charge(i) = 1.0
        sys%species_mass(i) = 1836.15
        call molecule_nullify_protons(sys%mol)
      else if (string_to_lowercase(sys%species_name(i)).eq.'positron') then
        sys%species_charge(i) = 1.0
        sys%species_mass(i) = 1.0
      else if (string_to_lowercase(sys%species_name(i)).eq.'muon') then
        sys%species_charge(i) = -1.0
        sys%species_mass(i) = 206.85
      else if (string_to_lowercase(sys%species_name(i)).eq.'nucleonp') then
        sys%species_charge(i) = 1.0
        sys%species_mass(i) = 1836.15
        call molecule_nullify_any(sys%mol)
      else if (string_to_lowercase(sys%species_name(i)).eq.'nucleonn') then
        sys%species_charge(i) = 0.0
        sys%species_mass(i) = 1838.68
        call molecule_nullify_any(sys%mol)
      else
        call file_error('Unknown species type')
      end if
      sys%moocc(1:sys%nmoo(i),i) = 1.0
      sys%moocc(sys%nmoo(i)+1:,i) = 0.0

      call get_basis(sys,sys%basis(i),basfile)

      if (sys%nmoo(i).gt.basis_nfun(sys%basis(i))) then
        call file_error('Negative number of virtual orbitals')
      end if
    end do
 1300 continue

    ! Interaction
    call set_default_interactions(sys)
    rewind(tape)
    read(tape,'(a)',end=1400) line
    do while(string_to_lowercase(line(1:13)).ne.'&interactions')
      read(tape,'(a)',end=1400) line
    end do
    read(tape,'(a)',end=1400) line
    do while(line(1:1).ne.'/')
      read(line,*) name1,name2
      i = index(line,'simple') 
      if (i.ne.0) then
        i = i + 7
        read(line(i:80),*) name3,val
        call assign_interaction(sys,name1,name2,'simple',name3,val)
      end if
      i = index(line,'exchange') 
      if (i.ne.0) then
        i = i + 9
        read(line(i:80),*) name3,val
        call assign_interaction(sys,name1,name2,'exchange',name3,val)
      end if
      i = index(line,'correlation') 
      if (i.ne.0) then
        i = i + 12
        read(line(i:80),*) name3,val
        call assign_interaction(sys,name1,name2,'correlation',name3,val)
      end if
      read(tape,'(a)',end=1400) line
    end do
 1400 continue

    ! Orbital type
    sys%orbital_type = 'cartesian' 
    rewind(tape)
    read(tape,'(a)',end=1500) line
    do while(string_to_lowercase(line(1:9)).ne.'&orbitals')
      read(tape,'(a)',end=1500) line
    end do
    read(tape,'(a)') sys%orbital_type
    sys%orbital_type = string_to_lowercase(sys%orbital_type)
 1500 continue

    return
99997 call file_error('Interaction block not found')
99998 call file_error('Basis set not found in input file')
99999 call file_error('Geometry block was not found in input file')
  end subroutine 

  integer function system_multiplicity(sys)
  implicit none
    type(nsystem) :: sys

    system_multiplicity = sys%nmoo(1) - sys%nmoo(2) + 1 
  end function

  real(8) function system_charge(sys)
  implicit none
    type(nsystem) :: sys

    integer :: iatom,is

    system_charge = 0.0
    do iatom=1,sys%mol%natom
      system_charge = system_charge + sys%mol%atom(iatom)%zeff
    end do
    do is=1,sys%ns
      system_charge = system_charge + &
      sys%species_charge(is)*float(sys%nmoo(is))
    end do
  end function

  subroutine get_basis(sys,basis,basfile)
  ! Load basis set according to input file
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys
    type(nbasis) :: basis
    character*(*) :: basfile

    integer :: iatom,iset,n

    basis%filename = basfile

    do iatom=1,sys%mol%natom
      iset = iatom
      n = sys%mol%atom(iatom)%atomic_number
      if ((n.eq.1).or.(n.eq.2)) then
        basis%set(iset)%firstrow = .true.
      else
        basis%set(iset)%firstrow = .false.
      end if
      call set_get(n,basis%set(iset),basfile) 
      basis%set(iset)%atom = iatom
      basis%nsets = basis%nsets + 1
    end do

    if (basis_nfun(basis).gt.nmaxnbas) &
      call file_error('Too many gaussian functions')
    return
99998 continue
    call file_error('Basis set not found')
  end subroutine

  subroutine get_ecps(sys,ecpfile)
  ! Load ECPs according to input file
  ! Roberto Flores-Moreno, 2020
  implicit none
    type(nsystem) :: sys
    character*(*) :: ecpfile

    integer :: iatom,iset,n

    sys%ecps%filename = ecpfile

    do iatom=1,sys%mol%natom
      iset = iatom
      n = sys%mol%atom(iatom)%atomic_number
      call set_ecp_get(n,sys%ecps%set(iset),ecpfile,sys%mol%atom(iatom)%zeff)
      sys%ecps%set(iset)%atom = iatom
      sys%ecps%nsets = sys%ecps%nsets + 1
    end do

    if (basis_nfun(sys%ecps).gt.nmaxnbas) &
      call file_error('Too many gaussian functions')
    return
99998 continue
    call file_error('ECPs not found')
  end subroutine

  subroutine set_default_interactions(sys)
  ! Set default type of exchange functional
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys

    integer :: is,js

    ! Default direct interactions based on particle types
    do is=1,sys%ns
      do js=is,sys%ns
        sys%interaction(is,js)%simple(:) = 0.0
        sys%interaction(is,js)%exchange(1) = 1.0  ! Fock exchange (4-ERIs)
        sys%interaction(is,js)%exchange(2) = 0.0
        sys%interaction(is,js)%correlation(1) = 0.0 
        if (((sys%species_name(is).eq.'electron').and.&
             (sys%species_name(js).eq.'electron')).or.&
            ((sys%species_name(is).eq.'proton').and.&
             (sys%species_name(js).eq.'proton'))) then
          sys%interaction(is,js)%simple(1) = 1.0    ! Coulomb repulsive (4-ERIs)
        else if (((sys%species_name(is).eq.'proton').and.&
                  (sys%species_name(js).eq.'electron')).or.&
                 ((sys%species_name(is).eq.'electron').and.&
                  (sys%species_name(js).eq.'proton'))) then
          sys%interaction(is,js)%simple(1) = -1.0  ! Coulomb attractive (4-ERIs)
        else if ((sys%species_name(is)(1:7).eq.'nucleon').and.&
                 (sys%species_name(js)(1:7).eq.'nucleon')) then
          sys%interaction(is,js)%simple(2) = 1.0   ! Yukawa
        end if
        sys%interaction(js,is) = sys%interaction(is,js)
      end do
    end do
  end subroutine

  subroutine assign_interaction(sys,names1,names2,namet,namei,val)
  ! Assign interaction for given pair of particle types
  ! Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: names1,names2,namet,namei
    real(8) :: val
    type(nsystem) :: sys
    
    integer :: is,js,n
    real(8) :: nval

    !rfm nval = max(0.0,min(val,1.0))
    nval = val

    do is=1,sys%ns
      do js=1,sys%ns
        if (((string_to_lowercase(names1).eq.&
              string_to_lowercase(sys%species_name(is))).and.&
             (string_to_lowercase(names2).eq.&
              string_to_lowercase(sys%species_name(js)))).or.&
            ((string_to_lowercase(names1).eq.&
              string_to_lowercase(sys%species_name(js))).and.&
             (string_to_lowercase(names2).eq.&
              string_to_lowercase(sys%species_name(is))))) then
          if (string_to_lowercase(namet).eq.'simple') then
            if (string_to_lowercase(namei).eq.'coulomb') then
              sys%interaction(is,js)%simple(1) = nval*sys%species_charge(is)*&
                                                     sys%species_charge(js)
            else if (string_to_lowercase(namei).eq.'none') then
              sys%interaction(is,js)%simple(:) = 0.0
            else if (string_to_lowercase(namei).eq.'yukawa') then
              sys%interaction(is,js)%simple(2) = nval
            else 
              print *, 'Unknown simple interaction type'
            end if
          else if (string_to_lowercase(namet).eq.'exchange') then
            n = xf_code(namei)
            if (n.le.0) then
              print *, 'Unknown exchange interaction type'
            else
              sys%interaction(is,js)%exchange(n) = nval
            end if
          else if (string_to_lowercase(namet).eq.'correlation') then
            n = cf_code(namei)
            if (n.le.0) then
              print *, 'Unknown correlation interaction type'
            else
              sys%interaction(is,js)%correlation(n) = nval
            end if
          end if
        end if
      end do
    end do
  end subroutine


  subroutine system_print_interactions(tape,sys)
  ! Print list of interactions
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: is,js
    real(8) :: w

    write(tape,*)
    do is=1,sys%ns
      do js=is,sys%ns
        write(tape,'(/,t2,"Interactions for pair: ",i4,2x,i4)') is,js
        write(tape,'(t2,"Simple: ")') 
        w = sys%interaction(is,js)%simple(1)
        if (w.ne.0.0) write(tape,'(t4,"Coulomb: ",f10.3)') w
        w = sys%interaction(is,js)%simple(2)
        if (w.ne.0.0) write(tape,'(t4,"Yukawa: ",f10.3)') w
        if (is.eq.js) then
          write(tape,'(t2,"Exchange: ")') 
          w = sys%interaction(is,js)%exchange(1)
          if (w.ne.0.0) write(tape,'(t4,"Fock: ",f10.3)') w
          w = sys%interaction(is,js)%exchange(2)
          if (w.ne.0.0) write(tape,'(t4,"Dirac: ",f10.3)') w
        end if
        write(tape,'(t2,"Correlation: ")') 
        w = sys%interaction(is,js)%correlation(1)
        if (w.ne.0.0) write(tape,'(t4,"VWN: ",f10.3)') w
      end do
    end do
    write(tape,*)
  end subroutine

  subroutine system_print(tape,sys)
  ! Print to tape system info
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: iatom,is
    real(8) :: charge

    ! Print geometry
    write(tape,'(t2,"Number of atoms: ",I6)') sys%mol%natom
    write(tape,'(t2,"Number of species types: ",I6)') sys%ns

    call file_string(tape,1,'Geometry in angstroms:',0)
    call molecule_print_geometry(tape,sys%mol)

    charge = 0
    do iatom=1,sys%mol%natom
      charge = charge + sys%mol%atom(iatom)%zeff
    end do
    do is=1,sys%ns
      charge = charge + sys%nmoo(is)*sys%species_charge(is)
    end do
    write(tape,'(t2,"Charge: ",F10.2)') charge
    write(tape,'(/,t2,"Species by type: ")') 
    do is=1,sys%ns
      write(tape,'(t4,"Type: ",i6,x,a20)') is,sys%species_name(is)
      write(tape,'(t6,"Number: ",i6)') sys%nmoo(is)
      write(tape,'(t6,"Mass: ",f20.10)') sys%species_mass(is)
      write(tape,'(t6,"Charge: ",f5.2)') sys%species_charge(is)
      write(tape,'(t6,"Number of basis functions: ",I6)') &
        basis_nfun(sys%basis(is))
      write(tape,'(t6,"Number of molecular orbitals: ",I6)') &
        basis_nfun(sys%basis(is))
      ! Print basis
      write(tape,'(//,t2,"Basis sets:")')
      call basis_print(tape,sys%basis(is))
      ! ECPs
      if (basis_nfun(sys%ecps).gt.0) then
        write(tape,'(//,t2,"ECPs:")')
        call basis_print(tape,sys%ecps)
      end if 
    end do
    
    write(tape,'(/,t2,"Interactions: ")') 
    call system_print_interactions(tape,sys)
  end subroutine

end module
