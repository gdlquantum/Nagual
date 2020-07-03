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
module nmgradient
! Evaluates forces
!
! R. Flores-Moreno, (PhD thesis, Cinvestav-IPN, 2006), equation (5.23)

  use nmparameter
  use nmtypes
  use nmfile
  use nmbasis
  use nmintegrals
  use nmshell
  use nmscf
  use nmadft
  use nmmatrix
  use nmxf
  use nmcf

  implicit none

    public :: gradient
    public :: gradient_print
    public :: gradient_pulay
    public :: gradient_kinetic
    public :: gradient_core
    public :: gradient_pulaym
    public :: gradient_kineticm
    public :: gradient_corem

    private

contains


  subroutine gradient(tape,sys,grad)
  ! Driver routine for gradient evaluation
  ! Roberto Flores-Moreno, 2008, 2010, 2018
  implicit none
    type(nsystem) :: sys
    integer :: tape
    real(8) :: grad(*)

    integer :: is,js
    logical :: dft

    integer :: allocs,dabas,dabasb,dbas
    real(8),allocatable :: c(:,:),p(:,:)
    real(8),allocatable :: xa(:),xb(:),yb(:),zb(:)

    dft = .false.
    if (xf_jacob(sys)+cf_jacob(sys).gt.0) dft = .true.
    grad(1:3*sys%mol%natom) = 0.0
 
    ! Nuclear repulsion gradients
    call getnregrad(sys%mol,grad)

    do is=1,sys%ns
      dbas = basis_nfun(sys%basis(is))
      dabas = basis_nfun(sys%auxis(is))
      allocate(c(dbas,dbas),p(dbas,dbas),xa(dabas),stat=allocs)
      if (allocs.gt.0) call file_error("gradient: allocation failed")

      ! Pulay forces
      call file_matrix_bas(c,dbas,is,'orbitals','read')
      call scf_build_pulay_matrix(sys,c,p,dbas,is)
      call gradient_pulay(sys%mol,sys%basis(is),p,grad,dbas,is)

      ! core and coulomb the full density matrix is used
      call file_matrix_bas(p,dbas,is,'density','read')
      call gradient_kinetic(sys,p,grad,dbas,is)
      call gradient_core(sys,p,grad,dbas,is)

      ! Coulomb/direct
      call file_vector_aux(xa,dabas,is,0,'x','read')
      do js=1,sys%ns
        call integrals_set_ppi_type(sys%interaction(is,js)%simple)

        dabasb = basis_nfun(sys%auxis(js))
        allocate(xb(dabasb),yb(dabasb),zb(dabasb),stat=allocs)
        if (allocs.gt.0) call file_error("gradient: allocation failed")

        call file_vector_aux(xb,dabasb,js,0,'x','read')
        yb(1:dabasb) = -0.5*xb(1:dabasb)
        if (dft.and.(js.eq.is)) then
          call file_vector_aux(zb,dabasb,js,0,'z','read')
          yb(1:dabasb) = yb(1:dabasb) - zb(1:dabasb)
        end if
        call gradient_dunlap(sys,xa,yb,grad,is,js)

        yb(1:dabasb) = xb(1:dabasb)
        if (dft.and.(js.eq.is)) then
          call file_vector_aux(zb,dabasb,js,0,'z','read')
          yb(1:dabasb) = yb(1:dabasb) + zb(1:dabasb)
        end if
        call gradient_aux_coulomb(sys,p,yb,grad,dbas,is,js)

        deallocate(xb,yb,zb,stat=allocs)
        if (allocs.gt.0) call file_error("gradient: deallocation failed")
      end do

      deallocate(c,p,xa,stat=allocs)
      if (allocs.gt.0) call file_error("gradient: deallocation failed")
    end do

    ! Exchange-correlation
    if (dft) call adft_xcegrad(sys,grad)
  end subroutine

  subroutine getnregrad(m,grad)
  ! Get nuclear repulsion energy gradients
  ! Roberto Flores-Moreno, 2009
  implicit none
    type(nmolecule) :: m
    real(8) :: grad(3,*)

    integer :: iatom,jatom
    real(8) :: ab,za,zb,ra(3),rb(3)

    do iatom=1,m%natom
      za = m%atom(iatom)%zeff
      ra = m%atom(iatom)%pos
      do jatom=iatom+1,m%natom
        zb = m%atom(jatom)%zeff
        rb = m%atom(jatom)%pos
        rb = rb - ra
        ab = sqrt(sum(rb(:3)**2))
        rb = za*zb*rb/ab**3
        grad(:3,iatom) = grad(:3,iatom) + rb(:3)
        grad(:3,jatom) = grad(:3,jatom) - rb(:3)
      end do
    end do
  end subroutine 

  subroutine gradient_add(la,lb,gblk,pblk,iblk,dblk,option)
  ! Apropiately pack gradients
  ! Roberto Flores-Moreno, Oct 2008
  implicit none
    character*(*) :: option
    integer :: dblk,la,lb
    real(8) :: gblk(3),pblk(dblk,dblk),iblk(dblk,dblk)

    integer :: a,ax,ay,az,b,bx,by,bz

    if (option.eq.'up') then
      do bx=0,lb
        do by=0,lb-bx
          bz = lb-bx-by
          b = raop(bx,by,bz)
          do ax=0,la
            do ay=0,la-ax
              az = la-ax-ay
              a = raop(ax,ay,az)
              gblk(1) = gblk(1) + 2.0*pblk(a,b)*iblk(raop(ax+1,ay,az),b)
              gblk(2) = gblk(2) + 2.0*pblk(a,b)*iblk(raop(ax,ay+1,az),b)
              gblk(3) = gblk(3) + 2.0*pblk(a,b)*iblk(raop(ax,ay,az+1),b)
            end do
          end do
        end do
      end do
    else if (option.eq.'down') then
      do bx=0,lb
        do by=0,lb-bx
          bz = lb-bx-by
          b = raop(bx,by,bz)
          do ax=0,la
            do ay=0,la-ax
              az = la-ax-ay
              a = raop(ax,ay,az)
              if (ax.gt.0) gblk(1) = gblk(1) - ax*                      &
     &                               pblk(a,b)*iblk(raop(ax-1,ay,az),b)
              if (ay.gt.0) gblk(2) = gblk(2) - ay*                      &
     &                               pblk(a,b)*iblk(raop(ax,ay-1,az),b)
              if (az.gt.0) gblk(3) = gblk(3) - az*                      &
    &                                pblk(a,b)*iblk(raop(ax,ay,az-1),b)
            end do
          end do
        end do
      end do
    end if
  end subroutine

  subroutine gradient_save(grad,gblk,atoml,atomr,atomc)
  ! Apropiately pack gradients
  ! Roberto Flores-Moreno, Oct 2008
  implicit none
    integer :: atoml,atomr,atomc
    real(8) :: grad(3,*),gblk(*)

    !!! Direct contributions
    grad(1:3,atoml) = grad(1:3,atoml) + gblk(1:3)

    !!! Traslational invariance contributions

    ! Two center integrals: B from A
    if (atomc.eq.0) then
      grad(1:3,atomr) = grad(1:3,atomr) - gblk(1:3)
    ! Three center integrals: C from A, C from B
    else
      grad(1:3,atomc) = grad(1:3,atomc) - gblk(1:3)
    end if
  end subroutine

  subroutine gradient_pulay(m,b,w,grad,dbas,is)
  ! Evaluates Pulay energy gradients
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: dbas,is
    real(8) :: w(dbas,dbas),grad(*)
    type(nmolecule) :: m
    type(nbasis) :: b

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,na,nb
    real(8) :: ra(3),rb(3)
    real(8) :: pblk(nmaxnco,nmaxnco),iblk(nmaxnco,nmaxnco),gblk(3)
    type(nshell) :: sa,sb,wsa

    ba = 0
    do iset=1,b%nsets
      ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,b%nsets
          rb(1:3) = m%atom(b%set(jset)%atom)%pos(1:3)
          gblk(1:3) = 0.0
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            pblk(1:na,1:nb) = w(ba+1:ba+na,bb+1:bb+nb)
            do i=1,na
              pblk(i,1:nb) = -pblk(i,1:nb)*sa%norm(i)
            end do
            ! Left up
            wsa%l = sa%l + 1
            call integrals_overlap(ra,rb,wsa,sb,iblk,nmaxnco,1)
            call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'up')
            ! Left down
            if (sa%l.gt.0) then
              wsa%l = sa%l - 1
              call integrals_overlap(ra,rb,wsa,sb,iblk,nmaxnco,0)
              call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'down')
            end if
            bb = bb + nb
          end do
          call gradient_save(grad,gblk,iatom,jatom,0)
        end do
       ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_kinetic(sys,p,grad,dbas,is)
  ! Evaluates kinetic energy gradients
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: dbas,is
    real(8) :: p(dbas,dbas),grad(*)
    type(nsystem) :: sys

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,na,nb
    real(8) :: ra(3),rb(3)
    real(8) :: pblk(nmaxnco,nmaxnco),iblk(nmaxnco,nmaxnco),gblk(3)
    type(nshell) :: sa,sb,wsa

    ba = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          gblk(1:3) = 0.0
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            pblk(1:na,1:nb) = p(ba+1:ba+na,bb+1:bb+nb)
            do i=1,na
              pblk(i,1:nb) = pblk(i,1:nb)*sa%norm(i)
            end do
            ! Left up
            wsa%l = sa%l + 1
            call integrals_kinetic(ra,rb,wsa,sb,iblk,1)
            call gradient_add(sa%l,sb%l,gblk,pblk,&
                iblk/sys%species_mass(is),nmaxnco,'up')
            ! Left down
            if (sa%l.gt.0) then
              wsa%l = sa%l - 1
              call integrals_kinetic(ra,rb,wsa,sb,iblk,0)
              call gradient_add(sa%l,sb%l,gblk,pblk,&
                               iblk/sys%species_mass(is),nmaxnco,'down')
            end if
            bb = bb + nb
          end do
          call gradient_save(grad,gblk,iatom,jatom,0)
        end do
        ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_core(sys,p,grad,dbas,is)
  ! Evaluates nuclear atraction energy gradients
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: dbas,is
    real(8) :: p(dbas,dbas),grad(*)
    type(nsystem) :: sys

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,katom,na,nb
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: pblk(nmaxnco,nmaxnco),iblk(nmaxnco,nmaxnco),gblk(3)
    type(nshell) :: sa,sb,wsa,wsb

    ba = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            wsb = sb 
            wsb%norm(:) = 1.0
            pblk(1:na,1:nb) = p(ba+1:ba+na,bb+1:bb+nb)
            do i=1,na
              do j=1,nb
                pblk(i,j) = pblk(i,j)*sa%norm(i)*sb%norm(j)
              end do
            end do
            do katom=1,sys%mol%natom
              if (sys%mol%atom(katom)%zeff.ne.0.0) then
                rc(1:3) = sys%mol%atom(katom)%pos(1:3)

                if (iatom.ne.katom) then
                  gblk(1:3) = 0.0
                  ! Up left
                  wsa%l = sa%l + 1
                  wsb%l = sb%l
                  call integrals_core(ra,rb,rc,wsa,wsb,iblk,1,0)
                  call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'up')
                  ! Down left
                  if (sa%l.gt.0) then
                    wsa%l = sa%l - 1
                    wsb%l = sb%l
                    call integrals_core(ra,rb,rc,wsa,wsb,iblk,0,0)
                    call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'down')
                  end if
                  gblk(1:3) = gblk(1:3)*sys%mol%atom(katom)%zeff*&
                                        sys%species_charge(is)
                  call gradient_save(grad,gblk,iatom,jatom,katom)
                end if

                if (jatom.ne.katom) then
                  gblk(1:3) = 0.0
                  pblk = transpose(pblk)
                  ! Up right
                  wsa%l = sa%l
                  wsb%l = sb%l + 1
                  call integrals_core(rb,ra,rc,wsb,wsa,iblk,1,0)
                  call gradient_add(sb%l,sa%l,gblk,pblk,iblk,nmaxnco,'up')
                  ! Down right
                  if (sb%l.gt.0) then
                    wsa%l = sa%l
                    wsb%l = sb%l - 1
                    call integrals_core(rb,ra,rc,wsb,wsa,iblk,0,0)
                    call gradient_add(sb%l,sa%l,gblk,pblk,iblk,nmaxnco,'down')
                  end if
                  pblk = transpose(pblk)
                  gblk(1:3) = gblk(1:3)*sys%mol%atom(katom)%zeff*&
                                        sys%species_charge(is)
                  call gradient_save(grad,gblk,jatom,iatom,katom)
                end if
              end if
            end do
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_print(tape,m,grad)
  ! Print out gradients
  ! Roberto Flores-Moreno, 2013
  implicit none
    type(nmolecule) :: m
    integer :: tape
    real(8) :: grad(3,*)

    integer :: i,iatom

    write(tape,'(/,t2,"Molecular gradients: ",/)') 
    do iatom=1,m%natom
      write(tape,'(t4,i6,x,3f15.5)') iatom,(grad(i,iatom),i=1,3)
    end do
  end subroutine 

  subroutine gradient_aux_coulomb(sys,p,x,grad,dbas,is,js)
  ! Evaluates fitted Coulomb energy gradients
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    integer :: dbas,is,js
    real(8) :: p(dbas,dbas),x(*),grad(*)
    type(nsystem) :: sys

    integer :: ba,bb,bc,bcc,i,iatom,iset,ishell,j,jatom,jset,jshell
    integer :: k,katom,kset,kshell,na,nb,nc
    logical :: skip
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: pblk(nmaxnco,nmaxnco),iblk(nmaxnco,nmaxnco),gblk(3)
    real(8) :: tblk(nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc,wsa,wsb

    ba = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = ((sa%l+1)*(sa%l+2))/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = ((sb%l+1)*(sb%l+2))/2
            wsb = sb 
            wsb%norm(:) = 1.0
            pblk(1:na,1:nb) = p(ba+1:ba+na,bb+1:bb+nb)
            do i=1,na
              do j=1,nb
                pblk(i,j) = pblk(i,j)*sa%norm(i)*sb%norm(j)
              end do
            end do
            bcc = 0
            do kset=1,sys%auxis(js)%nsets
              rc(1:3) = sys%mol%atom(sys%auxis(js)%set(kset)%atom)%pos(1:3)
              gblk(1:3) = 0.0
              ! Up left
              wsa%l = sa%l + 1
              wsb%l = sb%l
              iblk = 0.0
              bc = bcc
              do kshell=1,sys%auxis(js)%set(kset)%nshell
                sc = sys%auxis(js)%set(kset)%shell(kshell)
                nc = (sc%l+1)*(sc%l+2)/2
                call integrals_ppi_three(ra,rb,rc,wsa,wsb,sc,tblk,ntolint,&  
                                         skip,1,0)
                do k=1,nc
                  iblk(:,:) = iblk(:,:) + tblk(:,:,k)*x(bc+k)
                end do
                bc = bc + nc
              end do
              call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'up')
              ! Down left
              if (sa%l.gt.0) then
                wsa%l = sa%l - 1
                wsb%l = sb%l
                iblk = 0.0
                bc = bcc
                do kshell=1,sys%auxis(js)%set(kset)%nshell
                  sc = sys%auxis(js)%set(kset)%shell(kshell)
                  nc = (sc%l+1)*(sc%l+2)/2
                  call integrals_ppi_three(ra,rb,rc,wsa,wsb,sc,tblk,&
                                           ntolint,skip,0,0)
                  do k=1,nc
                    iblk(:,:) = iblk(:,:) + tblk(:,:,k)*x(bc+k)
                  end do
                  bc = bc + nc
                end do
                call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'down')
              end if
              call gradient_save(grad,gblk,iatom,jatom,katom)
              pblk = transpose(pblk)
              gblk(1:3) = 0.0
              ! Up right
              wsa%l = sa%l
              wsb%l = sb%l + 1
              iblk = 0.0
              bc = bcc
              do kshell=1,sys%auxis(js)%set(kset)%nshell
                sc = sys%auxis(js)%set(kset)%shell(kshell)
                nc = (sc%l+1)*(sc%l+2)/2
                call integrals_ppi_three(ra,rb,rc,wsa,wsb,sc,tblk,&
                                         ntolint,skip,0,1)
                do k=1,nc
                  iblk(:,:) = iblk(:,:) + tblk(:,:,k)*x(bc+k)
                end do
                bc = bc + nc
              end do
              iblk = transpose(iblk)
              call gradient_add(sb%l,sa%l,gblk,pblk,iblk,nmaxnco,'up')
              ! Down right
              if (sb%l.gt.0) then
                wsa%l = sa%l
                wsb%l = sb%l - 1
                iblk = 0.0
                bc = bcc
                do kshell=1,sys%auxis(js)%set(kset)%nshell
                  sc = sys%auxis(js)%set(kset)%shell(kshell)
                  nc = (sc%l+1)*(sc%l+2)/2
                  call integrals_ppi_three(ra,rb,rc,wsa,wsb,sc,tblk,&
                                           ntolint,skip,0,0)
                  do k=1,nc
                    iblk(:,:) = iblk(:,:) + tblk(:,:,k)*x(bc+k)
                  end do
                  bc = bc + nc
                end do
                iblk = transpose(iblk)
                call gradient_add(sb%l,sa%l,gblk,pblk,iblk,nmaxnco,'down')
              end if
              call gradient_save(grad,gblk,jatom,iatom,katom)
              pblk = transpose(pblk)
              bcc = bc
            end do
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_dunlap(sys,xl,xr,grad,is,js)
  ! Derivatives of Dunlap integrals
  ! Roberto Flores-Moreno, 2008, 2010, 2018
  implicit none
    integer :: is,js
    real(8) :: xl(*),xr(*),grad(*)
    type(nsystem) :: sys

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,na,nb
    real(8) :: factora
    real(8) :: gblk(3),ra(3),rb(3)
    real(8) :: pblk(nmaxnco,nmaxnco),iblk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,wsa

    ba = 0
    do iset=1,sys%auxis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%auxis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%auxis(is)%set(iset)%nshell
        sa = sys%auxis(is)%set(iset)%shell(ishell)
        na = ((sa%l+1)*(sa%l+2))/2
        wsa = sa
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%auxis(js)%nsets
          rb(1:3) = sys%mol%atom(sys%auxis(js)%set(jset)%atom)%pos(1:3)
          gblk(1:3) = 0.0 
          do jshell=1,sys%auxis(js)%set(jset)%nshell
            sb = sys%auxis(js)%set(jset)%shell(jshell)
            nb = ((sb%l+1)*(sb%l+2))/2
            do i=1,na
              factora =  xl(ba+i)*sa%norm(i)
              do j=1,nb
                pblk(i,j) = factora*xr(bb+j)
              end do
            end do
            wsa%l = sa%l + 1
            call integrals_ppi2(ra,rb,wsa,sb,iblk,1)
            call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'up')
            if (sa%l.gt.0) then
              wsa%l = sa%l - 1
              call integrals_ppi2(ra,rb,wsa,sb,iblk,0)
              call gradient_add(sa%l,sb%l,gblk,pblk,iblk,nmaxnco,'down')
            end if
            bb = bb + nb
          end do
          call gradient_save(grad,gblk,iatom,jatom,0)
        end do
        ba = ba + na
     end do
    end do
  end subroutine
 
  subroutine gradient_addm(la,lb,fa,fb,mat,iblk,norma,dbas,dblk,cc,option)
  ! Apropiately pack gradients for single matrix derivative
  ! Roberto Flores-Moreno, 2019
  implicit none
    character*(*) :: option
    integer :: cc,dbas,dblk,fa,fb,la,lb
    real(8) :: norma(*)
    real(8) :: mat(dbas,dbas),iblk(dblk,dblk)

    integer :: a,ax,ay,az,b,bx,by,bz
    real(8) :: factor

    if (option.eq.'up') then
      do bx=0,lb
        do by=0,lb-bx
          bz = lb-bx-by
          b = raop(bx,by,bz)
          do ax=0,la
            do ay=0,la-ax
              az = la-ax-ay
              a = raop(ax,ay,az)
              factor = 2.0*norma(a)
              if (cc.eq.1) then
                mat(fa+a,fb+b) = mat(fa+a,fb+b)+ factor*iblk(raop(ax+1,ay,az),b)
              else if (cc.eq.2) then
                mat(fa+a,fb+b) = mat(fa+a,fb+b)+ factor*iblk(raop(ax,ay+1,az),b)
              else if (cc.eq.3) then
                mat(fa+a,fb+b) = mat(fa+a,fb+b)+ factor*iblk(raop(ax,ay,az+1),b)
              end if
            end do
          end do
        end do
      end do
    else if (option.eq.'down') then
      do bx=0,lb
        do by=0,lb-bx
          bz = lb-bx-by
          b = raop(bx,by,bz)
          do ax=0,la
            do ay=0,la-ax
              az = la-ax-ay
              a = raop(ax,ay,az)
              if (cc.eq.1) then
                if (ax.gt.0) then
                  factor = ax*norma(a)
                  mat(fa+a,fb+b) = mat(fa+a,fb+b) - &
                                   factor*iblk(raop(ax-1,ay,az),b)
                end if
              else if (cc.eq.2) then
                if (ay.gt.0) then
                  factor = ay*norma(a)
                  mat(fa+a,fb+b) = mat(fa+a,fb+b) - &
                                   factor*iblk(raop(ax,ay-1,az),b)
                end if
              else if (cc.eq.3) then
                if (az.gt.0) then
                  factor = az*norma(a)
                  mat(fa+a,fb+b) = mat(fa+a,fb+b) - &
                                   factor*iblk(raop(ax,ay,az-1),b)
                end if
              end if
            end do
          end do
        end do
      end do
    end if
  end subroutine

  subroutine gradient_pulaym(m,b,s1,dbas,is,atom,cc)
  ! Evaluates Pulay energy gradients
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: atom,cc,dbas,is
    real(8) :: s1(dbas,dbas)
    type(nmolecule) :: m
    type(nbasis) :: b

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,na,nb
    real(8) :: ra(3),rb(3)
    real(8) :: iblk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,wsa

    s1(:,:) = 0.0

    ba = 0
    do iset=1,b%nsets
      ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,b%nsets
          rb(1:3) = m%atom(b%set(jset)%atom)%pos(1:3)
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if (((iatom.eq.atom).or.(jatom.eq.atom)).and.&
                (jatom.ne.iatom)) then
              ! Left up
              wsa%l = sa%l + 1
              call integrals_overlap(ra,rb,wsa,sb,iblk,nmaxnco,1)
              if (atom.eq.jatom) iblk = -iblk
              call gradient_addm(sa%l,sb%l,ba,bb,s1,iblk,sa%norm,&
                                 dbas,nmaxnco,cc,'up')
              ! Left down
              if (sa%l.gt.0) then
                wsa%l = sa%l - 1
                call integrals_overlap(ra,rb,wsa,sb,iblk,nmaxnco,0)
                if (atom.eq.jatom) iblk = -iblk
                call gradient_addm(sa%l,sb%l,ba,bb,s1,iblk,sa%norm,dbas,&
                                   nmaxnco,cc,'down')
              end if
            end if
            bb = bb + nb
          end do
        end do
       ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_kineticm(sys,t1,dbas,is,atom,cc)
  ! Evaluates kinetic energy gradients
  ! Roberto Flores-Moreno, 2008, 2009, 2018
  implicit none
    integer :: atom,cc,dbas,is
    real(8) :: t1(dbas,dbas)
    type(nsystem) :: sys

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,na,nb
    real(8) :: factor
    real(8) :: ra(3),rb(3)
    real(8) :: iblk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,wsa

    t1(1:dbas,1:dbas) = 0.0

    ba = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if (((iatom.eq.atom).or.(jatom.eq.atom)).and.&
                (jatom.ne.iatom)) then
              factor = 1.0/sys%species_mass(is)
              if (atom.eq.jatom) factor = -factor
              ! Left up
              wsa%l = sa%l + 1
              call integrals_kinetic(ra,rb,wsa,sb,iblk,1)
              iblk = iblk*factor
              call gradient_addm(sa%l,sb%l,ba,bb,t1,iblk,sa%norm,&
                                 dbas,nmaxnco,cc,'up')
              ! Left down
              if (sa%l.gt.0) then
                wsa%l = sa%l - 1
                call integrals_kinetic(ra,rb,wsa,sb,iblk,0)
                iblk = iblk*factor
                call gradient_addm(sa%l,sb%l,ba,bb,t1,&
                             iblk,sa%norm,dbas,nmaxnco,cc,'down')
              end if
            end if
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do
  end subroutine

  subroutine gradient_corem(sys,h1,dbas,is,atom,cc)
  ! Evaluates nuclear atraction matrix derivative
  ! Roberto Flores-Moreno, 2019
  implicit none
    integer :: atom,cc,dbas,is
    real(8) :: h1(dbas,dbas)
    type(nsystem) :: sys

    integer :: ba,bb,i,iatom,iset,ishell,ipt,j,jatom,jset,jshell,katom,na,nb
    real(8) :: factor
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: iblk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,wsa

    h1(1:dbas,1:dbas) = 0.0

    ba = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = ((sa%l+1)*(sa%l+2))/2
        wsa = sa 
        wsa%norm(:) = 1.0
        bb = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = ((sb%l+1)*(sb%l+2))/2
            do katom=1,sys%mol%natom
              if (((iatom.eq.atom).or.(katom.eq.atom)).and.&
                  (sys%mol%atom(katom)%zeff.ne.0.0)) then
                rc(1:3) = sys%mol%atom(katom)%pos(1:3)

                if (iatom.ne.katom) then
                  factor = sys%mol%atom(katom)%zeff*sys%species_charge(is)
                  if (atom.eq.katom) factor = -factor

                  ! Up 
                  wsa%l = sa%l + 1
                  call integrals_core(ra,rb,rc,wsa,sb,iblk,1,0)
                  iblk = iblk*factor
                  call gradient_addm(sa%l,sb%l,ba,bb,h1,&
                                     iblk,sa%norm,dbas,nmaxnco,cc,'up')
                  ! Down 
                  if (sa%l.gt.0) then
                    wsa%l = sa%l - 1
                    call integrals_core(ra,rb,rc,wsa,sb,iblk,0,0)
                    iblk = iblk*factor
                    call gradient_addm(sa%l,sb%l,ba,bb,h1,&
                                       iblk,sa%norm,dbas,nmaxnco,cc,'down')
                  end if
                end if
              end if
            end do
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do
    h1 = h1 + transpose(h1)
  end subroutine

end module


