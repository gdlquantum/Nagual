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
module nmmp
! Moller-Plesset perturbation theory driver
! Roberto Flores-Moreno, 2018

  use nmtypes
  use nmfile
  use nmbasis
  use nmaux
  use nmbecke

  implicit none

    public :: mp_driver

    private

      integer :: tape

contains

  subroutine mp_driver(tape,sys)
  ! Moller-Plesset perturbation theory driver
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: inp_tape,order
    namelist /mp/ order

    order = 0
    call file_open_ascii("Nagual.inp",inp_tape)
    read(inp_tape,nml=mp,end=1000)
1000 continue
    call file_close(inp_tape)
    if (order.le.1) return

    call aux_initialize(tape,sys)
    call aux_setup(tape,sys)
    call mp_2(tape,sys)
  end subroutine

  subroutine mp_2(tape,sys)
  ! MP2 energy
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: nq
    parameter (nq = 50)

    integer :: tape
    type(nsystem) :: sys

    real(8) :: mp2e
    real(8) :: rq(nq),wq(nq),buf(nq)

    ! Initialize 
    call file_string(tape,1,'=== MP ===',1)

    ! quadrature points
    call becke_radial_quadrature(rq,wq,nq)

    call mp_2_int(tape,sys,buf,rq,nq)
    mp2e = sum(buf(1:nq)*wq(1:nq))

    ! Add correlation energy to system energy
    sys%energy = sys%energy + mp2e

    write(tape,'(t2,"MP2 correlation energy: ",f20.6)') mp2e
    call file_flush(tape)
  end subroutine 

  subroutine mp_2_int(tape,sys,buf,rq,nq)
  ! Integration for MP2
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys
    integer :: nq,tape
    real(8) :: buf(*),rq(*)

    integer :: a,b,fa,i,is,iq,j,js,na,nmoo,nmoop,q
    real(8) :: aaw,e,s,t
    real(8) :: eig(nmaxnbas),eigp(nmaxnbas)

    integer :: allocs,dabas,dbas,dbasp,dorb,dorbp,dset
    logical,allocatable :: ldfaux(:),ldfbas(:)
    real(8),allocatable :: c(:,:),cp(:,:),eri(:,:),qeri(:,:),peri(:,:)
    real(8),allocatable :: g(:,:),x(:)
    real(8),allocatable :: lf(:,:),lfp(:,:)
   
    buf(1:nq) = 0.0
    aaw = 0.5
    dset = sys%mol%natom

    do js=1,sys%ns
      dbasp = basis_nfun(sys%basis(js))
      nmoop = sys%nmoo(js)
      eigp(1:dbasp) = sys%moe(1:dbasp,js)
      if (sys%orbital_type.ne.'cartesian') then
        dorbp = basis_nfun_spherical(sys%basis(js))
      else
        dorbp = dbasp
      end if
      allocate(cp(dbasp,dbasp),lfp(dorbp,nq),ldfaux(dset),ldfbas(dset),&
      stat=allocs)
      if (allocs.gt.0) call file_error('mp_2_int: allocation failed')

      do iq=1,nq
        do q=1,dorbp
          lfp(q,iq) = exp(-abs(eigp(q))*rq(iq))
        end do
      end do

      do j=1,nmoop
        ! Generate auxiliary basis for local density fitting
        call aux_ldf_pointer(sys,js,j,ldfaux,ldfbas,cp,dbasp)

        dabas = basis_restricted_nfun(sys%auxis(js),ldfaux)
        allocate(peri(dbasp,dabas),g(dabas,dabas),x(dabas),stat=allocs)
        if (allocs.gt.0) call file_error('mp_2_int: allocation failed')

        call file_matrix_bas(cp,dbasp,js,'orbitals','read')

        ! 3 index eris  ( p mu | k )
        call integrals_set_ppi_type(sys%interaction(js,js)%simple)
        call aux_mixed_eris(tape,cp(1,j),sys%mol,sys%basis(js),&
                       sys%auxis(js),ldfaux,peri,dabas,dbasp)

        ! 2 index eris (k | l) for weigthing 3 index eris
        call aux_pack_wg_matrix(js,sys%auxis(js),ldfaux,g,dabas)
        peri = matmul(peri,g)
        peri = matmul(peri,g)

        deallocate(g,x,stat=allocs)
        if (allocs.gt.0) call file_error('mp_2_int: deallocation failed')

        do is=1,js!,sys%ns
          dbas = basis_nfun(sys%basis(is))
          if (sys%orbital_type.ne.'cartesian') then
            dorb = basis_nfun_spherical(sys%basis(is))
          else
            dorb = dbas
          end if
          nmoo = sys%nmoo(is)
          eig(1:dbas) = sys%moe(1:dbas,is)
          allocate(qeri(dbas,dabas),eri(dbasp,dbas),lf(dorb,nq),&
               c(dbas,dbas),stat=allocs)
          if (allocs.gt.0) call file_error('mp_2_int: allocation failed')

          ! Laplace transformation factors
          do iq=1,nq
            do q=1,dorb
              lf(q,iq) = exp(-abs(eig(q))*rq(iq))
            end do
          end do

          call file_matrix_bas(c,dbas,is,'orbitals','read')

          call integrals_set_ppi_type(sys%interaction(is,js)%simple)
          do i=1,nmoo
            call aux_mixed_eris(tape,c(1,i),sys%mol,sys%basis(is),&
                           sys%auxis(js),ldfaux,qeri,dabas,dbas)
            eri = matmul(peri,transpose(qeri))
            eri = matmul(transpose(cp),eri)
            eri = matmul(eri,c)
            do a=nmoop+1,dorbp
              do b=nmoo+1,dorb
                if (is.eq.js) then
                  e = (eri(a,b) - eri(b,a))*aaw
                else
                  e = eri(a,b) 
                end if
                do iq=1,nq
                  t = e*e*lfp(j,iq)*lfp(a,iq)*lf(b,iq)*lf(i,iq)
                  buf(iq) = buf(iq) - t
                end do
              end do
            end do
          end do
          deallocate(qeri,eri,lf,c,stat=allocs)
          if (allocs.gt.0) call file_error('mp_2_int: deallocation failed')
        end do
        deallocate(peri,stat=allocs)
        if (allocs.gt.0) call file_error('mp_2_int: deallocation failed')
      end do
      deallocate(lfp,cp,ldfaux,ldfbas,stat=allocs)
      if (allocs.gt.0) call file_error('mp_2_int: deallocation failed')
    end do
  end subroutine

end module

