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
module nmg1pp
! Generalized one-particle propagator theory module
!
! Quasi-Particle G1PP
! Roberto Flores-Moreno, 2008, 2009
!
! R. Flores-Moreno, V. G. Zakrzewski and J. V. Ortiz, 
! J. Chem. Phys. 127, 134106 (2007).
!
! J. Romero, E. F. Posada, R. Flores-Moreno, A. Reyes,
! J. Chem. Phys. 137, 074105 (2012).
!
! M. A. Diaz-Tinoco, J. Romero, J. V. Ortiz, A. Reyes, R. Flores-Moreno,
! J. Chem. Phys. 138(19), 194108 (2013).
!
! Assumed ordering of spin-orbitals: p(alpha) h(alpha) p(beta) h(beta)
!

  use nmparameter
  use nmtypes
  use nmfile
  use nmsystem
  use nmtime
  use nmmolecule
  use nmintegrals
  use nmunits
  use nmmatrix
  use nmscf
  use nmbasis
  use nmset

  implicit none

    public :: g1pp_driver

    private

      integer :: tape

contains

  subroutine g1pp_driver(tape,sys,species_id,orbital_id,method,memory)
  ! Electron propagator theory driver
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    character*(*) :: method
    integer :: memory,orbital_id,species_id,tape
    type(nsystem) :: sys
    type(ntimer) :: tused

    real(8) :: pe,ps

    call file_header(tape,"G1PP calculation")

    call start_timer(tused,'Quasiparticle G1PP time:')
    call qpept(tape,sys,pe,ps,species_id,orbital_id,method,memory)
    call print_timer(tused,tape)
  end subroutine

  subroutine qpept(tape,sys,pe,ps,id,oid,method,memory)
  ! Quasi-Particle Electron Propagator Theory
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    integer :: nq
    parameter (nq = 15)

    character*(*) :: method   
    integer :: id,memory,oid,tape
    real(8) :: pe,ps
    type(nsystem) :: sys

    integer :: dbas,i,imo,is,rreq
    real(8) :: rq(nq),wq(nq),f2ph(nq),f2hp(nq)

    ! Initialize 
    call file_string(tape,1,'=== G1PP: Diagonal approximation ===',1)
    write(tape,'(t2,"Method: ",(a))') method
    write(tape,'(t2,"Species: ",i2)') id
    write(tape,'(t2,"Orbital: ",i4)') oid
    write(tape,'(t2,"Occupation: ",f10.2)') sys%moocc(oid,id)

    write(tape,'(t2,"Binding energies for species type ",i2)') id
    call file_string(tape,0,'==========================================',0)
    call file_string(tape,0,'  MO     Method       P.E.(eV)       P.S. ',0)
    call file_string(tape,0,'------------------------------------------',0)
    call file_flush(tape)

    ! KT
    pe = sys%moe(oid,id)
    ps = 1.0
    write(tape,5000) oid,'G1PP1',units_au_to_ev(pe),ps
    call file_flush(tape)
    if (method.eq.'d1') return

    call newton_raphson(sys,oid,pe,ps,id,selfe_2,memory)
    write(tape,5100) 'G1PP2',units_au_to_ev(pe),ps
    call file_flush(tape)
    if (method.eq.'d2') return 

    return
 5000 format(t2,i4,5x,a5,2x,2f13.3)
 5100 format(t11,a5,2x,2f13.3)
  end subroutine 

  subroutine newton_raphson(sys,p,pe,ps,id,selfe,memory)
  ! Newton-Raphson root search for QP EPT
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    integer :: id,memory,p
    real(8) :: pe,ps
    type(nsystem) :: sys
    external :: selfe

    integer :: it,topit
    real(8) :: corr,ds,pekt,ppe,s,tol

    ! settings
    topit = 20
    pekt = sys%moe(p,id)
    tol = 0.01/units_au_to_ev(1.0)

    pe = pekt
    do it=1,topit
      print *, pe
      ppe = pe
      call selfe(sys,id,p,ppe,s,ds,memory)
      print *, "S DS ",s,ds
      ps = 1.0/(1.0-ds)
      print *, "PS ",ps
      pe = ps*(pekt+s-ds*ppe)
      print *, "PEKT ",pekt
      corr = pe - ppe
      if (abs(corr).lt.tol) then
        corr = pe - pekt - s
        if (abs(corr).gt.tol) continue
        return
      end if
    end do
    write(tape,'(t2,"Newton-Raphson did not converge")')
  end subroutine 

  !!! Self-energy routines 
  subroutine selfe_2(sys,id,p,w,s,ds,memory)
  ! Second order self-energy for QP EPT
  ! 2008-2019, Roberto Flores-Moreno
  implicit none
    type(nsystem) :: sys
    integer :: id,memory,p
    real(8) :: w,s,ds(*)

    integer :: a,b,i,is,istep,j,llb,nmoo,nmoop,nstep,q,qq,ulb
    real(8) :: aaw,e,t,d

    integer :: allocs,dbas,dbasp,dorb,dorbp,dmo
    real(8),allocatable :: ii(:,:,:),c(:,:),cp(:,:),eig(:),eigp(:)

    ! initialize
    aaw = sqrt(0.5)
    s = 0.0
    ds(:2) = 0.0

    dbasp = basis_nfun(sys%basis(id))
    if (sys%orbital_type.ne.'cartesian') then
      dorbp = basis_nfun_spherical(sys%basis(id))
    else
      dorbp = dbasp
    end if
    allocate(cp(dbasp,dbasp),eigp(dbasp),stat=allocs)
    if (allocs.gt.0) call file_error('g1pp::selfe_2: allocation failed')
    nmoop = sys%nmoo(id)
    eigp(1:dbasp) = sys%moe(1:dbasp,id)
    call file_matrix_bas(cp,dbasp,id,'orbitals','read')

    do is=1,sys%ns
      dbas = basis_nfun(sys%basis(is))
      if (sys%orbital_type.ne.'cartesian') then
        dorb = basis_nfun_spherical(sys%basis(is))
      else
        dorb = dbas
      end if
      dmo = min(int(125000*memory/(dbas*dbasp)),dorb)
      print *, "SELFE_2, DMO/DBAS = ",dmo/float(dorb)
      if (dmo.lt.1) &
        call file_error('g1pp::selfe_2, memory is too low, aborting')
      allocate(ii(dbasp,dbas,dmo),c(dbas,dbas),eig(dbas),stat=allocs)
      if (allocs.gt.0) call file_error('g1pp::selfe_2: allocation failed')
      nmoo = sys%nmoo(is)
      eig(1:dbas) = sys%moe(1:dbas,is)
      call file_matrix_bas(c,dbas,is,'orbitals','read')

      nstep = dorb/dmo
      if (mod(dorb,dmo).ne.0) nstep = nstep + 1
      ulb = 0
      do istep=1,nstep
        llb = ulb + 1
        ulb = min(llb + dmo-1,dorb)
        call ii4batch(sys,p,id,is,ii,c,cp,dbas,dbasp,dorb,dorbp,llb,ulb)
        do q=llb,ulb
          qq = q - llb + 1
          ! 2ph
          if (q.le.nmoo) then
            do a=nmoop+1,dorbp
              do b=nmoo+1,dorb
                d = w + eig(q) - eigp(a) - eig(b)
                if (abs(d).gt.ntolnum) then ! Avoid zero in denominator
                  if (is.eq.id) then
                    e = (ii(a,b,qq) - ii(b,a,qq))*aaw
                  else
                    e = ii(a,b,qq)
                  end if
                  t = e*e/d
                  s = s + t
                  ds(1) = ds(1) - t/d
                end if
              end do
            end do
          ! 2hp
          else
            do i=1,nmoop
              do j=1,nmoo
                d = w + eig(q) - eigp(i) - eig(j)
                if (abs(d).gt.ntolnum) then ! Avoid zero in denominator
                  if (is.eq.id) then
                    e = (ii(i,j,qq)-ii(j,i,qq))*aaw
                  else
                    e = ii(i,j,qq)
                  end if
                  t = e*e/d
                  s = s + t
                  ds(1) = ds(1) - t/d
                end if
              end do
            end do
          end if
        end do
      end do
      deallocate(ii,c,eig,stat=allocs)
      if (allocs.gt.0) call file_error('g1pp::selfe_2: deallocation failed')
    end do
    deallocate(cp,eigp,stat=allocs)
    if (allocs.gt.0) call file_error('g1pp::selfe_2: deallocation failed')
  end subroutine

  subroutine ii4batch(sys,p,id,is,ii,c,cp,dbas,dbasp,dorb,dorbp,llb,ulb)
  ! Batch of four center ERIs in MO representation
  ! Roberto Flores Moreno, 2019
  implicit none
    type(nsystem) :: sys
    integer :: dbas,dbasp,dorb,dorbp,id,is,llb,p,ulb
    real(8) :: c(dbas,dbas),cp(dbasp,dbasp),ii(dbasp,dbas,*)

    integer :: db,fa,fb,fc,fd,iset,jset,kset,lset,ishell,jshell,kshell,lshell
    integer :: i,ib,j,jb,k,ibas,jbas,kbas,lbas,na,nb,nc,nd,q,r,t
    logical :: skip(2)
    real(8) :: ss
    real(8) :: ra(3),rb(3),rc(3),rd(3)
    real(8) :: cblk(nmaxnco,nmaxnco,nmaxnco)
    real(8) :: blk(nmaxnco,nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc,sd

    call integrals_set_ppi_type(sys%interaction(id,is)%simple)

    db = ulb-llb+1
    ii(:,:,1:db) = 0.0
    fa = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        fb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            fc = 0
            do kset=1,sys%basis(id)%nsets
              rc(1:3) = sys%mol%atom(sys%basis(id)%set(kset)%atom)%pos(1:3)
              do kshell=1,sys%basis(id)%set(kset)%nshell
                sc = sys%basis(id)%set(kset)%shell(kshell)
                nc = (sc%l+1)*(sc%l+2)/2
                ! ( p mu | nu sigma )
                cblk(1:na,1:nb,1:nc) = 0.0
                fd = 0
                do lset=1,sys%basis(id)%nsets
                  rd(1:3) = sys%mol%atom(sys%basis(id)%set(lset)%atom)%pos(1:3)
                  do lshell=1,sys%basis(id)%set(lset)%nshell
                    sd = sys%basis(id)%set(lset)%shell(lshell)
                    nd = (sd%l+1)*(sd%l+2)/2
                    call integrals_ppi(ra,rb,rc,rd,&
                                       sa,sb,sc,sd,blk,ntolint,skip,0,0,0)
                      if (skip(1)) go to 1000
                    if (.not.skip(2)) then
                      do ibas=1,na
                        do jbas=1,nb
                          do kbas=1,nc
                            ss = 0.0
                            do lbas=1,nd
                              ss= ss + cp(fd+lbas,p)*blk(ibas,jbas,kbas,lbas)
                            end do
                            cblk(ibas,jbas,kbas) = cblk(ibas,jbas,kbas) + ss
                          end do
                        end do
                      end do
                    end if
                    fd = fd + nd
                  end do
                end do
                ! ( p mu | nu t )
                do ib=llb,ulb
                  jb = ib - llb + 1
                  do ibas=1,na
                    ii(fc+1:fc+nc,fb+1:fb+nb,jb)=ii(fc+1:fc+nc,fb+1:fb+nb,jb)+&
                    c(fa+ibas,ib)*transpose(cblk(ibas,1:nb,1:nc))
                  end do
                end do
                fc = fc + nc
              end do
            end do
 1000       continue
            fb = fb + nb
          end do
        end do
        fa = fa + na
      end do
    end do
    ! ( p q | s t ) = < p s | q t >
    do ib=llb,ulb
      jb = ib - llb + 1
      ii(:,:dorb,jb) = matmul(ii(:,:,jb),c(:,:dorb))
    end do
    call matrix_transpose(cp,dbasp)
    do ib=llb,ulb
      jb = ib - llb + 1
      ii(1:dorbp,1:dorb,jb) = matmul(cp(:dorbp,:),ii(:,:dorb,jb))
    end do
    call matrix_transpose(cp,dbasp)
  end subroutine

end module

