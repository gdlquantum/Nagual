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
module nmfile
! File I/O management for Nagual

  use nmparameter
  use nmstring

  implicit none

    public :: file_open_ascii
    public :: file_open_binary
    public :: file_close
    public :: file_delete
    public :: file_flush
    public :: file_string
    public :: file_header
    public :: file_copyright
    public :: file_error
    public :: file_matrix_bas
    public :: file_matrix_bas_1
    public :: file_matrix_aux
    public :: file_matrix_diis
    public :: file_vector_aux
    public :: file_vector_eri4
    public :: file_grid

    private 

      character*(1) :: label(9)
      integer :: nfio,logt
      integer :: fstream(nmaxnfio)

      save

      data nfio /5/   ! Left out 5 and 6 for standard streams
      data label /"1","2","3","4","5","6","7","8","9"/

contains

  subroutine file_log
  ! Open LOG file for sequential access
  ! Roberto Flores Moreno, (Oct 2014)
  implicit none

    if (nfio.eq.5) then
      nfio = nfio + 1
      logt = nfio
      open(file="Nagual.log",&
           unit=logt,status='unknown',access='sequential')
    end if
  end subroutine

  subroutine file_open_ascii(filename,tape)
  ! Open file for sequential access
  ! Roberto Flores Moreno, 2008
  implicit none
    character*(*) :: filename
    integer :: tape

    call file_log

    nfio = nfio + 1
    if (nfio.gt.nmaxnfio) call file_error('too many files used this time')
    tape = nfio
    open(file=filename,unit=tape,status='unknown',access='sequential')
  end subroutine

  subroutine file_open_binary(filename,tape,rl)
  ! Open file for unformatted direct access (real data assumed)
  ! Roberto Flores Moreno, 2008
  implicit none
    character*(*) :: filename
    integer :: tape,rl

    call file_log

    nfio = nfio + 1
    if (nfio.gt.nmaxnfio) call file_error('too many files used this time')
    tape = nfio
    open(file=filename,unit=tape,status='unknown',access='direct',&
         form='unformatted',recl=8*rl)
  end subroutine

  subroutine file_close(tape)
  ! Close sequential access file
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape

    close(unit=tape,status='keep')
    ! If this was last open free its channel
    if (tape.eq.nfio) nfio = nfio - 1  
  end subroutine

  subroutine file_delete(tape)
  ! Delete a file
  ! Roberto Flores Moreno, May 2009
  implicit none
    integer :: tape

    close(unit=tape,status='delete')
    ! If this was last open free its channel
    if (tape.eq.nfio) nfio = nfio - 1  
  end subroutine

  subroutine file_string(tape,before,string,after)
  ! Print out a string
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: after,before,tape
    character*(*) :: string
    integer :: i

    do i=1,before
      write(tape,*)
    end do
    write(tape,'(T2,A)') string
    do i=1,after
      write(tape,*)
    end do
    
    call file_flush(tape)
  end subroutine 

  subroutine file_header(tape,string)
  ! Print out a header
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape
    character*(*) :: string

    call file_string(tape,2,'=== '//string//' ===',2)
  end subroutine

  subroutine file_copyright(tape)
  ! Print out copyright
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape

    call file_string(tape,1,&
      '=====================================================================',0)
    call file_string(tape,0,&
      '=                                                                   =',0)
    call file_string(tape,0,&
      '=                       <<< Nagual 1 >>>                            =',0)
    call file_string(tape,0,&
      '=                                                                   =',0)
    call file_string(tape,0,&
      '= R. Flores-Moreno, H.N. Gonzalez-Ramirez, J.F.H. Lew-Yee,          =',0)
    call file_string(tape,0,&
      '= J.M. del Campo, J.J. Villalobos-Castro, J.A. Flores-Ramos         =',0)
    call file_string(tape,0,&
      '= J.A. Guerrero-Cruz, B. Zuniga-Gutierrez,                          =',0)
    call file_string(tape,0,&
      '=                                                                   =',0)
    call file_string(tape,0,&
      '= E-mail: r.flores@academicos.udg.mx                                =',0)
    call file_string(tape,0,&
      '= Guadalajara Jal., Mexico 2020                                                    =',0)
    call file_string(tape,0,&
      '=                                                                   =',0)
    call file_string(tape,0,&
      '=====================================================================',1)
  end subroutine

  subroutine file_flush(tape)
  ! Flush file stream
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: tape

    call flush(tape)
  end subroutine

  subroutine file_error(string)
  ! Print out an error message and stop the program
  ! Roberto Flores Moreno, 2008
  implicit none
    character*(*) :: string

    call file_log

    call file_string(logt,1,string,1)
    call file_string(logt,1,'<<< Aborting >>>',1)
    do while (nfio.ge.6)
      call file_close(nfio)
    end do
    stop 'Internal stop'
  end subroutine 

  subroutine file_matrix_bas(m,dm,is,mname,option)
  ! Write/read SCF matrices from file
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    character*(*) :: mname,option
    integer :: dm,is
    real(8) :: m(dm,dm)

    integer :: i,j,k,recn,tape

    if (dm.le.0) call file_error('file_matrix_bas with null dimension')

    ! Open file
    call file_open_binary('Nagual.bas.'//label(is),tape,dm**2)

    ! Delete the file
    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    if (mname.eq.'overlap') then
      recn = 1
    else if (mname.eq.'overlaph') then
      recn = 2
    else if (mname.eq.'ortho') then
      recn = 3
    else if (mname.eq.'core') then
      recn = 4
    else if (mname.eq.'density') then
      recn = 5
    else if (mname.eq.'fock') then
      recn = 6
    else if (mname.eq.'orbitals') then
      recn = 7
    else if (mname.eq.'oorbitals') then
      recn = 8
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)
    end if

    ! Close file
    call file_close(tape)
  end subroutine 

  subroutine file_matrix_bas_1(m,dm,is,mname,option)
  ! Write/read perturbation matrices of SCF size from file
  ! Roberto Flores-Moreno, 2019
  implicit none
    character*(*) :: mname,option
    integer :: dm,is
    real(8) :: m(dm,dm)

    integer :: i,j,k,recn,tape

    if (dm.le.0) call file_error('file_matrix_bas with null dimension')

    ! Open file
    call file_open_binary('Nagual.bas1.'//label(is),tape,dm**2)

    ! Delete the file
    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    if (mname.eq.'overlap') then
      recn = 1
    else if (mname.eq.'core') then
      recn = 2
    else if (mname.eq.'density') then
      recn = 3
    else if (mname.eq.'fock') then
      recn = 4
    else if (mname.eq.'pulay') then
      recn = 5
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)
    end if

    ! Close file
    call file_close(tape)
  end subroutine 

  subroutine file_matrix_aux(m,dm,is,js,mname,option)
  ! Write/read auxiliary matrices to/from check point file
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    character*(*) :: mname,option
    integer :: dm,is,js
    real(8) :: m(*)

    integer :: i,iptr,ll,recn,tape,ul,rl

    ! Open file
    rl = dm
    if (rl.le.0) call file_error('file_matrix_aux with null dimension')
    call file_open_binary('Nagual.auxm.'//label(is)//'-'//label(js),tape,rl)

    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    ! Determine record number
    if (mname.eq.'g') then
      recn = 1
    else if (mname.eq.'ig') then
      recn = 2
    else if (mname.eq.'wg') then
      recn = 3
    else if (mname.eq.'a') then
      recn = 4
    else if (mname.eq.'f') then
      recn = 5
    else if (mname.eq.'r') then
      recn = 6
    else
      call file_error('file_matrix_aux: unknown matrix')
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) (m(i),i=1,rl)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) (m(i),i=1,rl)
    end if

    ! Close file
    call file_close(tape)
  end subroutine

  subroutine file_matrix_diis(m,dm,is,mname,option,recn,ndiis)
  ! Write/read fock and residual matrices from file for diis
  ! Felipe Lew-Yee, Jorge Martin del Campo-Ramirez, Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: mname,option
    integer :: dm,is
    real(8) :: m(dm,dm)

    integer :: i,j,recn,tape
    integer :: ndiis !number of fock matrices to use 

    if(ndiis.gt.0) then
      recn = MOD(recn,ndiis)
      if(recn.eq.0) then
        recn=ndiis
      end if
    end if

    if (dm.le.0) call file_error('file_matrix_diis with null dimension')

    ! Open file
    if(mname.eq.'r') then
      call file_open_binary('Nagual.diisr.'//label(is),tape,dm**2)
      else if(mname.eq.'f') then
      call file_open_binary('Nagual.diisf.'//label(is),tape,dm**2)
    end if

    ! Delete the file
    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) ((m(i,j),i=1,dm),j=1,dm)
    end if

    ! Close file
    call file_close(tape)
  end subroutine

  subroutine file_vector_aux(v,dm,is,ip,vname,option)
  ! Write/read auxiliary vectors to/from check point file
  ! Roberto Flores-Moreno, Oct 2008
  implicit none
    character*(*) :: vname,option
    integer :: dm,is,ip
    real(8) :: v(dm)

    integer :: i,j,recn,tape,rl

    rl = dm
    if (rl.le.0) call file_error('file_vector_aux with null dimension')

    ! Open file
    call file_open_binary('Nagual.auxv.'//label(is),tape,rl)

    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    ! Determine record number
    if (vname.eq.'x') then
      recn = 1*(ip+1)
    else if (vname.eq.'z') then
      recn = 2*(ip+1)
    else if (vname.eq.'b') then
      recn = 3*(ip + 1)
    end if

    if (option.eq.'write') then
      write(tape,rec=recn) (v(i),i=1,dm)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) (v(i),i=1,dm)
    end if

    call file_close(tape)
  end subroutine


  subroutine file_vector_eri4(m,dm,mname,option,recn,is)
  ! Write/read 3 center ERIs from file
  ! Felipe Lew-Yee, Jorge Martin del Campo-Ramirez, Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: mname,option
    integer :: dm
    real(8) :: m(dm)
    integer :: is
    integer :: i,recn,tape
    integer :: size_batch

    ! Open file
     call file_open_binary_int('Nagual.eri4info.'//label(is),tape,1)

     ! Delete the file
    if (option.eq.'delete') then
      call file_delete(tape)
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) dm
      size_batch = dm
    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) size_batch
    else if (option.eq.'delete') then
      size_batch = 1
    end if

    ! Close file
    call file_close(tape)

    ! Open file
    call file_open_binary('Nagual.eri4.'//label(is),tape,size_batch)

    ! Delete the file
    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    ! Write data
    if (option.eq.'write') then
      write(tape,rec=recn) (m(i),i=1,dm)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) (m(i),i=1,size_batch)
    end if

    ! Close file
    call file_close(tape)
  end subroutine

  subroutine file_open_binary_int(filename,tape,rl)
  ! Open file for unformatted direct access (int data assumed)
  ! Felipe Lew-Yee, Jorge Martin del Campo-Ramirez, Roberto Flores-Moreno, 2018
  implicit none
    character*(*) :: filename
    integer :: tape,rl

    call file_log

    nfio = nfio + 1
    if (nfio.gt.nmaxnfio) call file_error('too many files used this time')
    tape = nfio
    open(file=filename,unit=tape,status='unknown',access='direct',&
         form='unformatted',recl=4*rl)
  end subroutine

  subroutine file_grid(grid,dgrid,atom,option)
  ! Write/read auxiliary vectors to/from check point file
  ! Roberto Flores-Moreno, Oct 2008
  implicit none
    character*(*) :: option
    integer :: atom,dgrid
    real(8) :: grid(*)

    integer :: i,j,recn,tape,rl

    rl = dgrid
    if (rl.le.0) call file_error('file_grid with null dimension')

    ! Open file
    call file_open_binary('Nagual.grid',tape,rl)

    if (option.eq.'delete') then
      call file_delete(tape)
      return
    end if

    ! Determine record number
    recn = atom

    if (option.eq.'write') then
      write(tape,rec=recn) (grid(i),i=1,dgrid)

    ! Read data
    else if (option.eq.'read') then
      read(tape,rec=recn) (grid(i),i=1,dgrid)
    end if

    call file_close(tape)
  end subroutine

end module 
