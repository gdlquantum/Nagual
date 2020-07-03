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
module nmaux
! Handling of auxiliary functions

  use nmparameter
  use nmtypes
  use nmfile
  use nmbasis
  use nmintegrals
  use nmmatrix
  use nmsort

  implicit none

    public :: aux_initialize
    public :: aux_setup
    public :: aux_j_vector
    public :: aux_build_fock_matrix
    public :: aux_mixed_eris
    public :: aux_ldf_pointer
    public :: aux_pack_wg_matrix

contains

  subroutine aux_initialize(tape,sys)
  ! Setup RI
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: iset,ibas,ishell,is,na
    real(8) :: ra(3)
    real(8) :: blk(nmaxnco,nmaxnco)
    type(nshell) :: sa

    ! Generate auxiliary functions
    do is=1,sys%ns
      call integrals_set_ppi_type(sys%interaction(is,is)%simple)
      do iset=1,sys%basis(is)%nsets
        if ((sys%basis(is)%set(iset)%nshell.ne.0).and.&
            (sys%auxis(is)%set(iset)%nshell.eq.0)) then
          ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
          call auxgen(tape,sys%basis(is)%set(iset),&
                           sys%auxis(is)%set(iset))
          do ishell=1,sys%auxis(is)%set(iset)%nshell
            sa = sys%auxis(is)%set(iset)%shell(ishell)
            sa%norm = 1.0
            na = (sa%l+1)*(sa%l+2)/2
            call integrals_ppi_two(ra,ra,sa,sa,blk,0.0,0)
            do ibas=1,na
              sa%norm(ibas) = 1.0/sqrt(blk(ibas,ibas))
            end do
            sys%auxis(is)%set(iset)%shell(ishell) = sa
          end do
        end if
        sys%auxis(is)%nsets = sys%basis(is)%nsets
      end do
    end do
  end subroutine

  subroutine aux_setup(tape,sys)
  ! Setup RI
  ! Roberto Flores-Moreno, 2008, 2009, 2014, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: is

    integer :: allocs,dabas
    real(8),allocatable :: g(:,:),x(:)

    ! Matrices
    do is=1,sys%ns
      call integrals_set_ppi_type(sys%interaction(is,is)%simple)
      dabas = basis_nfun(sys%auxis(is))
      allocate(g(dabas,dabas),x(dabas),stat=allocs)
      if (allocs.gt.0) call file_error("aux_setup: allocation failed")
      call aux_build_g_matrix(tape,sys%mol,sys%auxis(is),g,dabas)
      call file_matrix_aux(g,dabas**2,is,is,'g','write')
      call matrix_svd_power(g,x,dabas,-1.0,1.0e-8)
      call file_matrix_aux(g,dabas**2,is,is,'ig','write')
      call file_matrix_aux(g,dabas**2,is,is,'g','read')
      call matrix_svd_power(g,x,dabas,-0.5,1.0e-8)
      call file_matrix_aux(g,dabas**2,is,is,'wg','write')
      deallocate(g,x,stat=allocs)
      if (allocs.gt.0) call file_error("aux_setup: deallocation failed")
    end do
  end subroutine

  subroutine auxgen(tape,set,setaux)
  ! Generate auxiliary functions for atom (even-tempered)
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: tape
    type(nset) :: set,setaux

    integer :: i,iblock,ishell,l,lmax,lbmax,n,naux,nauxis,nblock
    real(8) :: r,zmax,zmin,zbmax,zbmin
    type(nshell) :: s

    ! Max Z for basis
    lbmax = 0
    zbmax = 0.0
    do ishell=1,set%nshell
      s = set%shell(ishell)
      if (s%l.gt.lbmax) lbmax = s%l
      do i=1,s%k
        zbmax = max(zbmax,s%z(i))
      end do
    end do
    zbmin = zbmax
    do ishell=1,set%nshell
      s = set%shell(ishell)
      do i=1,s%k
        zbmin = min(zbmin,s%z(i))
      end do
    end do

    r = 4.0
    nblock = 3
    ! Give some flexibility to H and He
    if (set%firstrow) then
      r = 5.0
      nblock = 2
    end if

    ! How many z?
    n = int(log(zbmax/zbmin)/log(r) + 0.5)

    ! Larger z
    zmin = zbmin
    zmax= zmin*r**(n-1)
    zmax= 2.0*zmax  ! Remember we are trying to cover a product space

    s%k = 1
    s%c(1) = 1.0
    s%z(1) = zmax*r
    naux = 0
    do iblock=1,nblock
      if (iblock.eq.1) then
        nauxis = max(1,n/nblock + mod(n,nblock))
      else
        nauxis = max(1,n/nblock)
      end if
      do i=1,nauxis
        s%z(1) = s%z(1)/r
        if (i.eq.1) s%z(1) = s%z(1)*(r+1.0)/r
        do l=0,2*(iblock-1)
          naux = naux + 1
          if (naux.gt.nmaxshell) call file_error("nmaxshell")
          s%n = 0
          s%l = l
          setaux%shell(naux) = s
        end do
        if (i.eq.1) s%z(1) = s%z(1)*r/(r+1.0)
      end do
    end do
    setaux%nshell = naux
    setaux%atom = set%atom
  end subroutine

  subroutine aux_build_g_matrix(tape,m,b,g,dg)
  ! Build Dunlap (G) matrix
  ! Roberto Flores-Moreno, 2008, 2009
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: b
    integer :: dg,tape
    real(8) :: g(dg,dg)

    type(nshell) :: sa,sb
    integer :: ba,bb,iset,ibas,ishell,jset,jbas,jshell,na,nb
    real(8) :: ra(3),rb(3)
    real(8) :: blk(nmaxnco,nmaxnco)

    ba = 0
    do iset=1,b%nsets
      ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        bb = 0
        do jset=1,b%nsets
          rb(1:3) = m%atom(b%set(jset)%atom)%pos(1:3)
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if ((jset.gt.iset).or.                                    &
     &          ((jset.eq.iset).and.(jshell.ge.ishell))) then
              ! Put all of them explicitly calculated
              call integrals_ppi2(ra,rb,sa,sb,blk,0)
              do ibas=1,na
                do jbas=1,nb
                  g(ba+ibas,bb+jbas) = blk(ibas,jbas)
                end do
              end do
            end if
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do
    ! Now fa is equal to total number of auxiiary basis functions
    call matrix_symmetrize(g,ba,ba,'uplow')
  end subroutine

  subroutine aux_j_vector(p,j,m,basis,auxis,dbas,dabas)
  ! Build Coulomb vector J for VFCP
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis,auxis
    integer :: dabas,dbas
    real(8) :: j(dabas)
    real(8) :: p(dbas,dbas)

    integer :: ba,bb,bc,iset,ibas,ishell,jset,jbas,jshell,kset,kbas,    &
     & kshell,lamax,lmax,na,nb,nc
    logical :: skip
    real(8) :: factor
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: blk(nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc

    integer :: allocs,dleft,dright,dl
    real(8),allocatable :: h(:,:),hl(:,:),v(:,:,:)

    lmax = basis_lmax(basis)
    lamax = basis_lmax(auxis)
    dl = 2*lmax + lamax
    dleft = psdaop(2*lmax)
    dright = psdaop(lamax)
    allocate(h(dleft,dright),hl(dleft,dleft),v(dleft,dright,0:dl),&
             stat=allocs)
    if (allocs.gt.0) call file_error('aux_j_vector, allocation failed')

    j(1:dabas) = 0.0

    bc = 0
    do kset=1,auxis%nsets
      rc(1:3) = m%atom(auxis%set(kset)%atom)%pos(1:3)
      do kshell=1,auxis%set(kset)%nshell
        sc = auxis%set(kset)%shell(kshell)
        nc = (sc%l+1)*(sc%l+2)/2
        ba = 0
        do iset=1,basis%nsets
          ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
          do ishell=1,basis%set(iset)%nshell
            sa = basis%set(iset)%shell(ishell)
            na = (sa%l+1)*(sa%l+2)/2
            bb = 0
            do jset=1,basis%nsets
              rb(1:3) = m%atom(basis%set(jset)%atom)%pos(1:3)
              do jshell=1,basis%set(jset)%nshell
                sb = basis%set(jset)%shell(jshell)
                nb = (sb%l+1)*(sb%l+2)/2
               !call integrals_ppi_three(ra,rb,rc,sa,sb,sc,blk,ntolint,skip,0,0)
                call integrals_ppi3(ra,rb,rc,sa,sb,sc,blk,ntolint,skip,0,0,&
                            h,hl,v,dleft,dright,dl)
                do ibas=1,na
                  do jbas=1,nb
                    factor = sa%norm(ibas)*sb%norm(jbas)*p(ba+ibas,bb+jbas)
                    do kbas=1,nc
                      j(bc+kbas) = j(bc+kbas) +  factor*blk(ibas,jbas,kbas)
                    end do
                  end do
                end do
                bb = bb + nb
              end do
            end do
            ba = ba + na
          end do
        end do
        do kbas=1,nc
          j(bc+kbas) = j(bc+kbas)*sc%norm(kbas)
        end do
        bc = bc + nc
      end do
    end do

    deallocate(h,hl,v,stat=allocs)
    if (allocs.gt.0) call file_error('aux_j_vector, deallocation failed')

    ! Normalize
  end subroutine

  subroutine aux_build_fock_matrix(tape,m,basis,auxis,f,x,dbas,dabas)
  ! Build Fock matrix using VFCP
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis,auxis
    integer :: dabas,dbas,is,js,tape
    real(8) :: x(dabas)
    real(8) :: f(dbas,dbas)

    integer :: ba,bb,bc,iset,ibas,ishell,jset,jbas,jshell,kset,kbas,    &
       kshell,lamax,lmax,na,nb,nc
    logical :: skip
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: blk(nmaxnco,nmaxnco,nmaxnco),sblk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc

    integer :: allocs,dleft,dright,dl
    real(8),allocatable :: h(:,:),hl(:,:),v(:,:,:),y(:)

    lmax = basis_lmax(basis)
    lamax = basis_lmax(auxis)
    dl = 2*lmax + lamax
    dleft = psdaop(2*lmax)
    dright = psdaop(lamax)
    allocate(h(dleft,dright),hl(dleft,dleft),v(dleft,dright,0:dl),&
             y(dabas),stat=allocs)
    if (allocs.gt.0) call file_error('aux_fock_matrix, allocation failed')

    bc = 0
    do kset=1,auxis%nsets
      rc(1:3) = m%atom(auxis%set(kset)%atom)%pos(1:3)
      do kshell=1,auxis%set(kset)%nshell
        sc = auxis%set(kset)%shell(kshell)
        nc = (sc%l+1)*(sc%l+2)/2
        do kbas=1,nc
          y(bc+kbas) = x(bc+kbas)*sc%norm(kbas)
        end do
        bc = bc + nc
      end do
    end do

    ba = 0
    do iset=1,basis%nsets
      ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
      do ishell=1,basis%set(iset)%nshell
        sa = basis%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        bb = 0
        do jset=1,basis%nsets
          rb(1:3) = m%atom(basis%set(jset)%atom)%pos(1:3)
          do jshell=1,basis%set(jset)%nshell
            sb = basis%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            sblk(1:na,1:nb) = 0.0
            bc = 0
            do kset=1,auxis%nsets
              rc(1:3) = m%atom(auxis%set(kset)%atom)%pos(1:3)
              do kshell=1,auxis%set(kset)%nshell
                sc = auxis%set(kset)%shell(kshell)
                nc = (sc%l+1)*(sc%l+2)/2
                !call integrals_ppi_three(ra,rb,rc,sa,sb,sc,blk,ntolint,skip,0,0)
                call integrals_ppi3(ra,rb,rc,sa,sb,sc,blk,ntolint,skip,0,0,&
                            h,hl,v,dleft,dright,dl)
                do kbas=1,nc
                  sblk(1:na,1:nb) = sblk(1:na,1:nb) + &
                  blk(1:na,1:nb,kbas)*y(bc+kbas)
                end do
                bc = bc + nc
              end do
            end do
            do jbas=1,nb
              do ibas=1,na
                f(ba+ibas,bb+jbas) = f(ba+ibas,bb+jbas) + sblk(ibas,jbas)*&
                sa%norm(ibas)*sb%norm(jbas)
              end do
            end do
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do

    deallocate(h,hl,v,y,stat=allocs)
    if (allocs.gt.0) call file_error('aux_fock_matrix, deallocation failed')
  end subroutine

  subroutine aux_mixed_eris(tape,cp,m,b,baux,ptr,rectm,dabas,dbas)
  ! Store three center ERIs 
  ! Robert Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: b,baux
    integer :: dabas,dbas,tape
    logical :: ptr(*)
    real(8) :: cp(*)
    real(8) :: rectm(dbas,dabas)

    integer :: atoma,atomb,ba,bb,bc,iset,ibas,ishell,jset,jbas,&
               jshell,kset,kbas,kshell,na,nb,nc
    logical :: skip
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: cblk(nmaxnco,nmaxnco),blk(nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc

    bc = 0
    do kset=1,baux%nsets
      if (ptr(kset)) then
        rc(1:3) = m%atom(baux%set(kset)%atom)%pos(1:3)
        do kshell=1,baux%set(kset)%nshell
          sc = baux%set(kset)%shell(kshell)
          nc = (sc%l+1)*(sc%l+2)/2
          ba = 0
          do iset=1,b%nsets
            ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
            do ishell=1,b%set(iset)%nshell
              sa = b%set(iset)%shell(ishell)
              na = (sa%l+1)*(sa%l+2)/2
              cblk(1:na,1:nc) = 0.0
              bb = 0
              do jset=1,b%nsets
                rb(1:3) = m%atom(b%set(jset)%atom)%pos(1:3)
                do jshell=1,b%set(jset)%nshell
                  sb = b%set(jset)%shell(jshell)
                  nb = (sb%l+1)*(sb%l+2)/2
                  call integrals_ppi_three(ra,rb,rc,sa,sb,sc,blk,ntolint,skip,0,0)
                  do jbas=1,nb
                    cblk(1:na,1:nc) = cblk(1:na,1:nc) + &
                    blk(1:na,jbas,1:nc)*cp(bb+jbas)
                  end do
                  bb = bb + nb
                end do
              end do
              rectm(ba+1:ba+na,bc+1:bc+nc) = cblk(1:na,1:nc) 
              ba = ba + na
            end do
          end do
          bc = bc + nc
        end do
      end if
    end do
  end subroutine

  subroutine aux_ldf_pointer(sys,id,oid,ptr,ptrb,c,dbas)
  ! Address auxiliary basis for local density fitting
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: dbas,id,oid
    logical :: ptr(*),ptrb(*)
    real(8) :: c(dbas,dbas)
    type(nsystem) :: sys

    integer :: dset,fa,iset,ishell,na
    integer :: setptr(nmaxset)
    real(8) :: s
    real(8) :: setw(nmaxset)
    type(nshell) :: sa

    ! Load orthogonalized MO coefficients
    call file_matrix_bas(c,dbas,id,'oorbitals','read')

    dset = sys%basis(id)%nsets

    fa = 0
    do iset=1,dset
      na = 0
      do ishell=1,sys%basis(id)%set(iset)%nshell
        sa = sys%basis(id)%set(iset)%shell(ishell)
        na = na + ((sa%l+1)*(sa%l+2))/2
      end do
      setw(iset) = sum(c(fa+1:fa+na,oid)**2)
      setptr(iset) = iset
      fa = fa + na
    end do

    call sort_vector(setw,setptr,dset)

    ptr(1:dset) = .false.
    ptrb(1:dset) = .false.
    s = 1.0
    do iset=dset,1,-1
      if (s.gt.1.0e-4) ptr(setptr(iset)) = .true.
      if (s.gt.1.0e-8) ptrb(setptr(iset)) = .true.
      s = s - setw(iset)
    end do
  end subroutine

  subroutine aux_pack_wg_matrix(is,b,ptr,g,dg)
  ! Pack LDF G matrix
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nbasis) :: b
    integer :: dg,is,tape
    logical :: ptr(*)
    real(8) :: g(dg,dg)

    type(nshell) :: sa,sb
    integer :: ba,bb,iset,ishell,jset,na,nb
    integer :: llset(nmaxset),ulset(nmaxset)
    
    integer :: allocs,dabas
    real(8),allocatable :: big(:,:),x(:)

    dabas = basis_nfun(b)
    if (dabas.eq.dg) then
      call file_matrix_aux(g,dg**2,is,is,'wg','read')
    else

      ba = 0
      do iset=1,b%nsets
        llset(iset) = ba + 1
        na = 0
        do ishell=1,b%set(iset)%nshell
          sa = b%set(iset)%shell(ishell)
          na = na + ((sa%l+1)*(sa%l+2))/2
        end do
        ba = ba + na
        ulset(iset) = ba
      end do

      allocate(big(dabas,dabas),x(dg),stat=allocs)
      if (allocs.gt.0) call file_error("aux_pack_g_matrix: allocation failed")

      call file_matrix_aux(big,dabas**2,is,is,'g','read')

      ba = 0
      do iset=1,b%nsets
        if (ptr(iset)) then
          na = ulset(iset) - llset(iset) + 1
          bb = ba
          do jset=iset,b%nsets
            if (ptr(jset)) then
              nb = ulset(jset) - llset(jset) + 1
              g(ba+1:ba+na,bb+1:bb+nb) = big(llset(iset):ulset(iset),&
                                             llset(jset):ulset(jset))
              bb = bb + nb
            end if
          end do
          ba = ba + na
        end if
      end do

      call matrix_symmetrize(g,dg,dg,'uplow')
      call matrix_svd_power(g,x,dg,-0.5,1.0e-8)

      deallocate(big,x,stat=allocs)
      if (allocs.gt.0) call file_error("aux_pack_g_matrix: deallocation failed")
    end if
  end subroutine

end module


