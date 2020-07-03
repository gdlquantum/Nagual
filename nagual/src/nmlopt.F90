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
module nmlopt
! Local optimizer

  use nmparameter
  use nmtypes
  use nmfile
  use nmmolecule
  use nmenergy
  use nmgradient
  use nmgui
  use nmxf
  use nmcf
  use nmgrid
  use nmelement
 
  implicit none

    public :: lopt

    private

contains

  subroutine lopt(ncycles,sys,task)
  ! Driver routine for the geometry optimization.
  ! Roberto Flores-Moreno, 2014, 2018
  implicit none
    integer :: ncycles
    type(ntask) :: task
    type(nsystem) :: sys

    integer :: i,iatom,nit,tape
    logical :: optimized
    real(8) :: energy,ng,gmax
    real(8) :: de,fscala,norm,step,penergy
    type(nmolecule) :: bkp

    integer :: allocs,datom
    real(8),allocatable :: grad(:,:),ngrad(:,:),oldgrad(:,:)

    call file_open_ascii("Nagual.opt",tape)
    call file_copyright(tape)
    call file_header(tape,'Local molecular geometry optimization')

    ! initialization 
    datom = sys%mol%natom
    allocate(grad(3,datom),ngrad(3,datom),oldgrad(3,datom),stat=allocs)
    if (allocs.gt.0) call file_error('opt_local: allocation failed')

    !call energy_single_point(tape,task,sys,energy,.false.)
    !call gradient(tape,sys,grad)
    !call test_gradient(tape,task,sys,grad,ngrad)
    !goto 3333
    fscala = 2.0
    nit = 0
    optimized = .false.
    do while (.not.optimized)
      nit = nit + 1

      call energy_single_point(tape,task,sys,energy,.false.)

      call gradient(tape,sys,grad)
      gmax = maxval(abs(grad(1:3,1:sys%mol%natom)))
      if (abs(gmax).lt.1.0e-3) optimized = .true.

      write (tape,'(t2,"------------------------------")')   
      write (tape,'(t2,"Local optimization cycle: ",i4)') nit 
      write (tape,'(t2,"Current energy: ",f20.9)') energy
      write (tape,'(t2,"Current geometry:")')
      call molecule_print_geometry(tape,sys%mol)

      write (tape,'(t2,"Current gradients:")')
      do iatom=1,sys%mol%natom
        write(tape,'(t2,a2,x,3f20.5)') sys%mol%atom(iatom)%symbol,&
         (grad(i,iatom),i=1,3)
      end do

      call gui_optimization_step(sys%mol,grad,energy,nit)

      norm = sqrt(sum(grad(1:3,1:sys%mol%natom)*grad(1:3,1:sys%mol%natom)))
      grad(1:3,1:sys%mol%natom) = grad(1:3,1:sys%mol%natom)/norm 

      step = 0.1
      penergy = energy
      de = 0.1
      do while (de.ge.0.0) 
        do iatom=1,sys%mol%natom
          sys%mol%atom(iatom)%pos(:) = sys%mol%atom(iatom)%pos(:) - &
            grad(:,iatom)*step
        end do
        call energy_single_point(tape,task,sys,energy,.false.)
        de = energy - penergy
        if (de.lt.0.0) exit
        do iatom=1,sys%mol%natom
          sys%mol%atom(iatom)%pos(:) = sys%mol%atom(iatom)%pos(:) + &
            grad(:,iatom)*step
        end do
        step = step/fscala
      end do

      ng = de/step
      write (tape,'(t2,"Step: ",f20.9)') step
      write (tape,'(t2,"Energy change: ",f20.9)') de
      write (tape,'(t2,"Numerical directional derivative: ",f20.9)') ng

      call file_flush(tape)
      if (nit.ge.ncycles) exit
    end do

    if (optimized) then
      write (tape,'(/,t2,"Geometry is OPTIMIZED")')
    else
      write (tape,'(/,t2,"Geometry is NOT optimized")')
    end if

!3333 continue
    deallocate(grad,ngrad,oldgrad,stat=allocs)
    if (allocs.gt.0) call file_error('opt_local: deallocation failed')

    call file_close(tape)
  end subroutine 

  subroutine test_gradient(tape,task,sys,agrad,ngrad)
  ! Testing gradients
  ! Roberto Flores-Moreno, 2019
  implicit none
    integer :: tape
    real(8) :: agrad(3,*),ngrad(3,*)
    type(ntask) :: task
    type(nsystem) :: sys

    integer :: cc,iatom,ih
    real(8) :: ep,em,shift

    ! Numerical
    shift = 1.0e-3
    write(tape,'(/,t1,"Testing gradients:")')
    write(tape,'(t1,"====================================")')
    ih = 0
    do iatom=1,sys%mol%natom
      write(tape,'(t1,"Atom ",i4)') iatom
      do cc=1,3
        ih = ih + 1
        sys%mol%atom(iatom)%pos(cc) = sys%mol%atom(iatom)%pos(cc) + shift
        task%scf%guess = 'core'
        call energy_single_point(tape,task,sys,sys%energy,.false.)
        ep = sys%energy
        sys%mol%atom(iatom)%pos(cc) = sys%mol%atom(iatom)%pos(cc) - 2.0*shift
        task%scf%guess = 'core'
        call energy_single_point(tape,task,sys,sys%energy,.false.)
        em = sys%energy
        sys%mol%atom(iatom)%pos(cc) = sys%mol%atom(iatom)%pos(cc) + shift
        ngrad(cc,iatom) = (ep - em)/(2.0*shift)
        write(tape,'(t1,i1,x,3f20.5)') cc,agrad(cc,iatom),ngrad(cc,iatom),&
                                       agrad(cc,iatom)-ngrad(cc,iatom)
      end do
    end do
    write(tape,'(t1,"====================================")')
  end subroutine

end module 

