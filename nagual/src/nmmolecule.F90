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
module nmmolecule
! Molecule 
  
  use nmparameter
  use nmtypes
  use nmfile
  use nmstring
  use nmatom
  use nmmatrix
  use nmelement
  use nmvector
  use nmunits

  implicit none

    public :: molecule_copy
    public :: molecule_print_geometry
    public :: molecule_read_cartesian
    public :: molecule_read_zmatrix
    public :: molecule_nullify_protons
    public :: molecule_nullify_any

    private

contains

  subroutine molecule_copy(src,dst)
  ! Copy object information
  ! Roberto Flores-Moreno (Oct 2014)
  implicit none
    type(nmolecule) :: dst,src

    integer :: iatom

    dst%pg = src%pg
    dst%natom = src%natom
    do iatom=1,src%natom
      call atom_copy(src%atom(iatom),dst%atom(iatom))
    end do
  end subroutine

  subroutine molecule_print_geometry(tape,m)
  ! Print geometry to a file
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape
    type(nmolecule) :: m

    integer :: iatom

    do iatom=1,m%natom
      call atom_print_geometry(tape,m%atom(iatom))
    end do
  end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Dimension reporters

  subroutine molecule_read_cartesian(m,tape,flag,inbohr)
  ! Load molecular cartesian coordinates
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    character :: flag
    integer :: tape
    logical :: inbohr

    character*(80) :: line
    character*(2) :: symbol
    integer :: i,natom

    natom = 0
    read(tape,'(a)') line
    do while(line(1:1).ne.flag)
      natom = natom + 1
      read(line,*) symbol,(m%atom(natom)%pos(i),i=1,3)
      call element_symbol_to_atomic_number(symbol,m%atom(natom)%atomic_number)
      m%atom(natom)%symbol = element_get_symbol(m%atom(natom)%atomic_number)
      m%atom(natom)%zeff = real(m%atom(natom)%atomic_number)
      ! Angstrom to Bohr
      if (.not.inbohr) then
        do i=1,3
          m%atom(natom)%pos(i) = units_angstrom_to_bohr(m%atom(natom)%pos(i))
        end do
      end if
      read(tape,'(a)') line
    end do
    m%natom = natom
  end subroutine

  subroutine molecule_read_zmatrix(m,tape,flag,inbohr)
  ! Load molecular Z-Matrix coordinates
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    character :: flag
    integer :: tape
    logical :: inbohr

    character*(80) :: line,str
    character*(2) :: symbol
    character*(10) :: atomlab(nmaxatom)
    character*(10) :: varllab(3*nmaxatom),reflab(3,nmaxatom),varlab(3,nmaxatom)
    integer :: aref,bref,c0,c9,dref,i,iatom,ivar,i2,jatom,natom,nref,nvar
    integer :: refval(3,nmaxatom)
    logical :: found,zmat
    real(8) :: aux(3),axis(3),shift(3)
    real(8) :: varlval(3*nmaxatom),varval(3,nmaxatom)

    ! reading
    zmat = .true.
    natom = 0
    nvar = 0
    read(tape,'(t1,a80)') line
    do while(line(1:1).ne.flag)
      if (line(1:1).eq.'#') then
        zmat = .false.
      else if (zmat) then
        natom = natom + 1
        nref = min(3,natom-1)
        if (nref.eq.0) then
          read(line,*) atomlab(natom)
        else if (nref.eq.1) then
          read(line,*) atomlab(natom),reflab(1,natom),varlab(1,natom)
        else
          read(line,*) atomlab(natom),(reflab(i,natom),varlab(i,natom),i=1,nref)
        end if
      else
        nvar = nvar + 1
        read(line,*) varllab(nvar),varlval(nvar)
      end if
      read(tape,'(t1,a80)') line
    end do
    m%natom = natom

    ! variable substitution
    c0 = ichar('0')
    c9 = ichar('9')
    do iatom=1,m%natom
      i2 = ichar(atomlab(iatom)(2:2))
      if ((i2.ge.c0).and.(i2.le.c9)) then 
        symbol(1:1) = atomlab(iatom)(1:1)
        symbol(2:2) = ' '
      else
        symbol(1:2) = atomlab(iatom)(1:2)
      end if
      call element_symbol_to_atomic_number(symbol,m%atom(iatom)%atomic_number)
      m%atom(iatom)%symbol = element_get_symbol(m%atom(iatom)%atomic_number)
      m%atom(iatom)%zeff = real(m%atom(iatom)%atomic_number)
      nref = min(3,iatom-1)
      do i=1,nref
        ! Look for reference values
        jatom = 1
        do while(jatom.lt.iatom)
          if (string_to_lowercase(reflab(i,iatom)).eq.&
              string_to_lowercase(atomlab(jatom))) exit
          jatom = jatom + 1
        end do
        if (jatom.ge.iatom) call file_error('MOLECULE: Z matrix ref. not found')
        refval(i,iatom) = jatom
        ! Look for variable values
        i2 = ichar(varlab(i,iatom)(1:1))
        if ((i2.ge.c0).and.(i2.le.c9)) then 
           read(varlab(i,iatom),*) varval(i,iatom)
        else
          found = .false.
          do ivar=1,nvar
            if ((string_to_lowercase(varlab(i,iatom)).eq.&
                 string_to_lowercase(varllab(ivar))).or.&
                (string_to_lowercase(varlab(i,iatom)(2:10)).eq.&
                 string_to_lowercase(varllab(ivar)(1:9)))) then
              varval(i,iatom) = varlval(ivar) 
              if (string_to_lowercase(varlab(i,iatom)(1:1)).eq.'-') then
                varval(i,iatom) = -varlval(ivar) 
              end if
              found = .true.
            end if
          end do
          if (.not.found) call file_error('MOLECULE: Z matrix var. not found')
        end if
        if ((i.eq.1).and.(.not.inbohr)) then
          varval(i,iatom) = units_angstrom_to_bohr(varval(i,iatom))
        end if
      end do
      ! Get cartesians
      m%atom(iatom)%pos(1:3) = 0.0
      if (iatom.gt.1) then
        bref = refval(1,iatom)
        m%atom(iatom)%pos(1:3) = m%atom(bref)%pos(1:3) 
        if (iatom.gt.2) then
          aref = refval(2,iatom)
          shift(1:3) = m%atom(aref)%pos(1:3) - m%atom(bref)%pos(1:3)
          call vector_normalize(shift)
        else
          call vector_constructor(shift,0.0,0.0,1.0)
        end if
        shift(1:3) = shift(1:3)*varval(1,iatom)
        m%atom(iatom)%pos(1:3) = m%atom(iatom)%pos(1:3) + shift(1:3)
        if (iatom.gt.2) then
          shift(1:3) = m%atom(aref)%pos(1:3) - m%atom(bref)%pos(1:3)
          if (iatom.eq.3) then
            call vector_constructor(aux,0.0,1.0,0.0)
            call vector_cross(aux,shift,axis)
         else
            dref = refval(3,iatom)
            aux(1:3) = m%atom(dref)%pos(1:3) - m%atom(bref)%pos(1:3)
            call vector_cross(shift,aux,axis)
          end if
          shift(1:3) = m%atom(iatom)%pos(1:3) - m%atom(bref)%pos(1:3)
          if (vector_colineal(shift,axis)) call vector_ortho(shift,axis)
          call vector_rotate(shift,axis,varval(2,iatom))
          m%atom(iatom)%pos(1:3) = m%atom(bref)%pos(1:3) + shift(1:3)
          if (iatom.gt.3) then
            axis(1:3) = m%atom(aref)%pos(1:3) - m%atom(bref)%pos(1:3)
            shift(1:3) = m%atom(iatom)%pos(1:3) - m%atom(bref)%pos(1:3) 
            call vector_rotate(shift,axis,varval(3,iatom))
            m%atom(iatom)%pos(1:3) = m%atom(bref)%pos(1:3) + shift(1:3)
          end if
        end if
      end if
    end do
  end subroutine

  subroutine molecule_nullify_protons(m)
  ! Remove classical point charges corresponding to protons
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m

    integer :: iatom

    do iatom=1,m%natom
      if (m%atom(iatom)%atomic_number.eq.1) then
        m%atom(iatom)%zeff = 0.0
      end if
    end do
  end subroutine

  subroutine molecule_nullify_any(m)
  ! Remove classical point charges corresponding to any atom
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m

    integer :: iatom

    do iatom=1,m%natom
      m%atom(iatom)%zeff = 0.0
    end do
  end subroutine

end module
