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
! B. Zuniga-Gutierrez, Nagual 1, Guadalajara, Mexico (2020)
!
!###################################################
!   Nagual: Multicomponent many body calculations.
!   Copyright (C) 2006-2020 Nagual developers.
!
!   Contact: r.flores@academicos.udg.mx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nmadft
! ADFT 

  use nmtypes
  use nmfile
  use nmstring
  use nmbecke
  use nmlebedev
  use nmbasis
  use nmset
  use nmshell
  use nmxf
  use nmcf
  use nmvector
  use nmtime
  use nmgrid


  implicit none

    public :: adft_vxc
    public :: adft_xcegrad

    private 

      real(8) :: minrho        ! Smallest density processed in grid (DFT)
      parameter (minrho = 1.0e-10)

contains

  subroutine adft_vxc(tape,sys,e,dft)
  ! EXchange--Correlation potential and energy
  ! Roberto Flores-Moreno, 2009, 2010, 2018
  implicit none
    integer :: tape
    logical :: dft
    real(8) :: e
    type(nsystem) :: sys

    integer :: i,iap,iatom,igp,irp,is,js
    type(ntimer) :: t1,t2,t3,t4

    integer :: allocs,dabas,ds
    real(8),allocatable :: g(:,:),z(:)

    e = 0.0
    if (.not.dft) return

    ds = sys%ns

    ! Exchange
    if (xf_jacob(sys).gt.0) then
      do is=1,ds
        call adft_vx(tape,sys,e,is)
      end do
    end if

    ! Correlation
    if (cf_jacob(sys).gt.0) then
      do is=1,ds
        do js=is+1,ds
          call adft_vc(tape,sys,e,is,js)
        end do
      end do
    end if

    ! Apply G^{-1}
    do is=1,ds
      dabas = basis_nfun(sys%auxis(is))
      allocate(g(dabas,dabas),z(dabas),stat=allocs)
      if (allocs.gt.0) call file_error("adft_vxc: Allocation failed")
      call file_vector_aux(z,dabas,is,0,'z','read')
      call file_matrix_aux(g,dabas**2,is,is,'ig','read')
      z = matmul(g,z)
      call file_vector_aux(z,dabas,is,0,'z','write')
      deallocate(g,z,stat=allocs)
      if (allocs.gt.0) call file_error("adft_vxc: Deallocation failed")
    end do
  end subroutine

  subroutine adft_vx(tape,sys,e,is)
  ! EXchange potential and energy
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,tape
    real(8) :: e
    type(nsystem) :: sys

    integer :: i,iap,iatom,iset,igp,irp,jset,ll,ul
    logical :: skip(nmaxset)
    real(8) :: drho,rho,v,ev,t

    integer :: allocs,dabas
    real(8),allocatable :: x(:),z(:)
    real(8),allocatable :: abas(:),gbuf(:,:)

    dabas = basis_nfun(sys%auxis(is))
    allocate(abas(dabas),z(dabas),x(dabas),gbuf(4,nmaxgpa),stat=allocs)
    if (allocs.gt.0) call file_error("adft_vx: Allocation failed")

    z = 0.0
    call file_vector_aux(x,dabas,is,0,'x','read')

    do iset=1,sys%auxis(is)%nsets
      iatom = iset
      call file_grid(gbuf,4*nmaxgpa,iset,'read')
      do igp=1,ngpc(iset)
        call evaluate_basis(sys%mol,sys%auxis(is),gbuf(1,igp),abas,dabas,skip)
        rho = sum(x(1:dabas)*abas(1:dabas))
        if (rho.gt.minrho) then
          call xf_energy(sys,rho,ev,is)
          e = e + gbuf(4,igp)*ev
          call xf_potential(sys,rho,v,is)
          v = gbuf(4,igp)*v
          z(1:dabas) = z(1:dabas)  + v*abas(1:dabas)
        end if
      end do
    end do

    call file_vector_aux(z,dabas,is,0,'z','write')

    deallocate(abas,z,x,gbuf,stat=allocs)
    if (allocs.gt.0) call file_error("adft_vx: Deallocation failed")
  end subroutine

  subroutine adft_vc(tape,sys,e,is,js)
  ! Correlation potential and energy
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,js,tape
    real(8) :: e
    type(nsystem) :: sys

    integer :: i,iap,iatom,igp,irp,iset,jset,ll,ul
    logical :: skipa(nmaxset),skipb(nmaxset)
    real(8) :: rhoa,rhob,vca,vcb,ev

    integer :: allocs,dabasa,dabasb
    real(8),allocatable :: xa(:),za(:),abasa(:),gbuf(:,:)
    real(8),allocatable :: xb(:),zb(:),abasb(:)

    dabasa = basis_nfun(sys%auxis(is))
    dabasb = basis_nfun(sys%auxis(js))
    allocate(abasa(dabasa),abasb(dabasb),za(dabasa),zb(dabasb),& 
             xa(dabasa),xb(dabasb),gbuf(4,nmaxgpa),stat=allocs)
    if (allocs.gt.0) call file_error("adft_vc: Allocation failed")

    call file_vector_aux(xa,dabasa,is,0,'x','read')
    call file_vector_aux(xb,dabasb,js,0,'x','read')
    call file_vector_aux(za,dabasa,is,0,'z','read')
    call file_vector_aux(zb,dabasb,js,0,'z','read')

    do iset=1,sys%auxis(is)%nsets
      iatom = iset
      call file_grid(gbuf,4*nmaxgpa,iset,'read')
      do igp=1,ngpc(iset)
        call evaluate_basis(sys%mol,sys%auxis(is),gbuf(1,igp),abasa,dabasa,skipa)
        call evaluate_basis(sys%mol,sys%auxis(js),gbuf(1,igp),abasb,dabasb,skipb)
        rhoa = sum(xa(1:dabasa)*abasa(1:dabasa))
        rhob = sum(xb(1:dabasb)*abasb(1:dabasb))
        if ((rhoa.gt.0.0).and.(rhob.gt.0.0).and.&
            (max(rhoa,rhob).gt.minrho)) then
          call cf_energy(sys,rhoa,rhob,ev,is,js)
          e = e + gbuf(4,igp)*ev
          call cf_potential(sys,rhoa,rhob,vca,vcb,is,js)
          vca = gbuf(4,igp)*vca
          vcb = gbuf(4,igp)*vcb
          za(1:dabasa) = za(1:dabasa)  + vca*abasa(1:dabasa)
          zb(1:dabasb) = zb(1:dabasb)  + vcb*abasb(1:dabasb)
        end if
      end do
    end do

    call file_vector_aux(za,dabasa,is,0,'z','write')
    call file_vector_aux(zb,dabasb,js,0,'z','write')

    deallocate(abasa,abasb,za,zb,xa,xb,gbuf,stat=allocs)
    if (allocs.gt.0) call file_error("adft_vc: Deallocation failed")
  end subroutine


  subroutine adft_xcegrad(sys,gradient)
  ! EXchange--Correlation energy gradient
  ! Roberto Flores-Moreno, 2018
  implicit none
    real(8) :: gradient(3,*)
    type(nsystem) :: sys

    integer :: is,js
    integer :: allocs,dabasa,dabasb
    real(8),allocatable :: xa(:),xb(:)

    ! Exchange
    if (xf_jacob(sys).gt.0) then
      do is=1,sys%ns
        dabasa = basis_nfun(sys%auxis(is))
        allocate(xa(dabasa),stat=allocs)
        if (allocs.gt.0) call file_error("adft_xcegrad: Allocation failed")
        call file_vector_aux(xa,dabasa,is,0,'x','read')
        call adft_xegrad(sys,gradient,is,xa,xa)
        deallocate(xa,stat=allocs)
        if (allocs.gt.0) call file_error("adft_xcegrad: Deallocation failed")
      end do
    end if

    ! Correlation
    if (cf_jacob(sys).gt.0) then
      do is=1,sys%ns
        dabasa = basis_nfun(sys%auxis(is))
        allocate(xa(dabasa),stat=allocs)
        if (allocs.gt.0) call file_error("adft_xcegrad: Allocation failed")
        call file_vector_aux(xa,dabasa,is,0,'x','read')
        do js=is+1,sys%ns
          dabasb = basis_nfun(sys%auxis(js))
          allocate(xb(dabasb),stat=allocs)
          if (allocs.gt.0) call file_error("adft_xcegrad: Allocation failed")
          call file_vector_aux(xb,dabasb,js,0,'x','read')
          call adft_cegrad(sys,gradient,is,js,xa,xb,xa,xb)
          deallocate(xb,stat=allocs)
          if (allocs.gt.0) call file_error("adft_xcegrad: Deallocation failed")
        end do
        deallocate(xa,stat=allocs)
        if (allocs.gt.0) call file_error("adft_xcegrad: Deallocation failed")
      end do
    end if
  end subroutine

  subroutine adft_xegrad(sys,gradient,is,xrho,xsum)
  ! EXchange energy gradients
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is
    real(8) :: gradient(3,*),xrho(*),xsum(*)
    type(nsystem) :: sys

    integer :: i,iap,iatom,iset,igp,irp,jatom,jset,ll,ul
    logical :: skip(nmaxset)
    real(8) :: rho,v,ev,t

    integer :: allocs,dabas
    real(8),allocatable :: abas(:),abasd1(:,:),gbuf(:,:)

    dabas = basis_nfun(sys%auxis(is))
    allocate(abas(dabas),abasd1(3,dabas),gbuf(4,nmaxgpa),stat=allocs)
    if (allocs.gt.0) call file_error("adft_vx_grad: allocation failed")

    do iset=1,sys%auxis(is)%nsets
      iatom = iset
      call file_grid(gbuf,4*nmaxgpa,iset,'read')
      do igp=1,ngpc(iset)
        call evaluate_basis_grad(sys%mol,sys%auxis(is),gbuf(1,igp),abas,&
                                 abasd1,dabas,skip)
        rho = sum(xrho(1:dabas)*abas(1:dabas))
        if (rho.gt.minrho) then
          call xf_potential(sys,rho,v,is)
          v = gbuf(4,igp)*v
          do jset=1,sys%auxis(is)%nsets
            if (.not.skip(jset)) then
              jatom = jset
              ll = sys%auxis(is)%set(jset)%ll
              ul = sys%auxis(is)%set(jset)%ul
              do i=1,3
                gradient(i,jatom) = gradient(i,jatom) + &
                v*sum(abasd1(i,ll:ul)*xsum(ll:ul))
              end do
            end if
          end do
        end if
      end do
    end do

    deallocate(abas,abasd1,gbuf,stat=allocs)
    if (allocs.gt.0) call file_error("adft_vx_grad: deallocation failed")
  end subroutine

  subroutine adft_cegrad(sys,gradient,is,js,xarho,xbrho,xasum,xbsum)
  ! Correlation energy gradient
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: is,js
    real(8) :: gradient(3,*),xarho(*),xbrho(*),xasum(*),xbsum(*)
    type(nsystem) :: sys

    integer :: i,iap,iatom,igp,irp,iset,jatom,jset,ll,ul
    logical :: skipa(nmaxset),skipb(nmaxset)
    real(8) :: rhoa,rhob,vca,vcb,ev

    integer :: allocs,dabasa,dabasb
    real(8),allocatable :: abasa(:),abasad1(:,:)
    real(8),allocatable :: abasb(:),abasbd1(:,:)
    real(8),allocatable :: gbuf(:,:)

    dabasa = basis_nfun(sys%auxis(is))
    dabasb = basis_nfun(sys%auxis(js))
    allocate(abasa(dabasa),abasad1(3,dabasa),&
             abasb(dabasb),abasbd1(3,dabasb),&
             gbuf(4,nmaxgpa),stat=allocs)
    if (allocs.gt.0) call file_error("adft_cegrad: allocation failed")

    do iset=1,sys%auxis(is)%nsets
      iatom = iset
      call file_grid(gbuf,4*nmaxgpa,iset,'read')
      do igp=1,ngpc(iset)
        call evaluate_basis_grad(sys%mol,sys%auxis(is),gbuf(1,igp),&
                                      abasa,abasad1,dabasa,skipa)
        call evaluate_basis_grad(sys%mol,sys%auxis(js),gbuf(1,igp),&
                                      abasb,abasbd1,dabasb,skipb)
        rhoa = sum(xarho(1:dabasa)*abasa(1:dabasa))
        rhob = sum(xbrho(1:dabasb)*abasb(1:dabasb))
        if (max(rhoa,rhob).gt.minrho) then
          call cf_potential(sys,rhoa,rhob,vca,vcb,is,js)
          vca = gbuf(4,igp)*vca
          do jset=1,sys%auxis(is)%nsets
            if (.not.skipa(jset)) then
              jatom = jset
              ll = sys%auxis(is)%set(jset)%ll
              ul = sys%auxis(is)%set(jset)%ul
              do i=1,3
                gradient(i,jatom) = gradient(i,jatom) + &
                vca*sum(abasad1(i,ll:ul)*xasum(ll:ul))
              end do
            end if
          end do
          vcb = gbuf(4,igp)*vcb
          do jset=1,sys%auxis(js)%nsets
            if (.not.skipb(jset)) then
              jatom = jset
              ll = sys%auxis(js)%set(jset)%ll
              ul = sys%auxis(js)%set(jset)%ul
              do i=1,3
                gradient(i,jatom) = gradient(i,jatom) + &
                vcb*sum(abasbd1(i,ll:ul)*xbsum(ll:ul))
              end do
            end if
          end do
        end if
      end do
    end do

    deallocate(abasa,abasb,abasad1,abasbd1,gbuf,stat=allocs)
    if (allocs.gt.0) call file_error("adft_cegrad: deallocation failed")
  end subroutine

  subroutine evaluate_basis(m,basis,r,abas,dabas,skip)
  ! Evaluates auxiliary functions
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    integer :: dabas
    logical :: skip(*)
    real(8) :: r(3),abas(dabas)
    type(nmolecule) :: m
    type(nbasis) :: basis

    integer :: ba,i,iatom,iset,irp,ishell,l,lmax,lx,ly,lz,na
    real(8) :: expf,r2,z,zp
    real(8) :: ra(3),p(3,0:nmaxl)
    type(nshell) :: sa

    abas(1:dabas) = 0.0

    p(:3,0) = 1.0
    do iset=1,basis%nsets
      iatom = iset
      ra(1:3) = r(1:3) - m%atom(iatom)%pos(1:3)
      r2 = sum(ra*ra)
      if (r2.lt.basis%set(iset)%radius) then
        skip(iset) = .false.
        ba = basis%set(iset)%ll - 1
        lmax = set_lmax(basis%set(iset))
        do l=1,lmax
          p(1:3,l) = p(1:3,l-1)*ra(1:3)
        end do
        zp = 0.0
        do ishell=1,basis%set(iset)%nshell
          sa = basis%set(iset)%shell(ishell)
          l = sa%l
          na = shell_nco(l)
          z = sa%z(1)
          if (z.ne.zp) then
            expf = exp(-z*r2)
            zp = z
          end if
          do lx=l,0,-1
            do ly=l-lx,0,-1
              lz = l - lx - ly
              i = raop(lx,ly,lz) 
              abas(ba+i) = expf*sa%norm(i)*p(1,lx)*p(2,ly)*p(3,lz)
            end do
          end do
          ba = ba + na
        end do
        skip(iset) = .false.
      else
        skip(iset) = .true.
      end if
    end do
  end subroutine

  subroutine evaluate_basis_grad(m,basis,r,abas,abasd1,dabas,skip)
  ! Evaluates auxiliary functions
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: dabas
    logical :: skip(*)
    real(8) :: r(3),abas(dabas),abasd1(3,dabas)
    type(nmolecule) :: m
    type(nbasis) :: basis

    integer :: ba,bai,i,iatom,iset,irp,ishell,l,lmax,lx,ly,lz,na
    real(8) :: factor,expf,r2,z,z2,zp
    real(8) :: ra(3),p(3,0:nmaxl+1)
    type(nshell) :: sa

    abasd1(3,1:dabas) = 0.0

    p(:3,0) = 1.0
    do iset=1,basis%nsets
      iatom = iset
      ra(1:3) = r(1:3) - m%atom(iatom)%pos(1:3)
      r2 = sum(ra*ra)
      if (r2.lt.basis%set(iset)%radius) then
        skip(iset) = .false.
        ba = basis%set(iset)%ll - 1
        lmax = set_lmax(basis%set(iset)) + 1  
        do l=1,lmax
          p(1:3,l) = p(1:3,l-1)*ra(1:3)
        end do
        zp = 0.0
        do ishell=1,basis%set(iset)%nshell
          sa = basis%set(iset)%shell(ishell)
          l = sa%l
          na = shell_nco(l)
          z = sa%z(1)
          if (z.ne.zp) then
            expf = exp(-z*r2)
            zp = z
          end if
          z2 = 2.0*z
          do lx=l,0,-1
            do ly=l-lx,0,-1
              lz = l - lx - ly
              i = raop(lx,ly,lz) 
              bai = ba + i
              abasd1(1,bai) = z2*p(1,lx+1)
              abasd1(2,bai) = z2*p(2,ly+1)
              abasd1(3,bai) = z2*p(3,lz+1)
              if (lx.gt.0) abasd1(1,bai) = abasd1(1,bai) - float(lx)*p(1,lx-1)
              if (ly.gt.0) abasd1(2,bai) = abasd1(2,bai) - float(ly)*p(2,ly-1)
              if (lz.gt.0) abasd1(3,bai) = abasd1(3,bai) - float(lz)*p(3,lz-1)
              abasd1(1,bai) = abasd1(1,bai)*p(2,ly)*p(3,lz)
              abasd1(2,bai) = abasd1(2,bai)*p(1,lx)*p(3,lz)
              abasd1(3,bai) = abasd1(3,bai)*p(2,ly)*p(1,lx)
              factor = expf*sa%norm(i)
              abasd1(1:3,bai) = abasd1(1:3,bai)*factor
              abas(bai) = factor*p(1,lx)*p(2,ly)*p(3,lz)
            end do
          end do
          ba = ba + na
        end do
        skip(iset) = .false.
      else
        skip(iset) = .true.
      end if
    end do
  end subroutine

end module


