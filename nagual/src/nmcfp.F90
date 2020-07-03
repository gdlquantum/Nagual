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
module nmcfp

  use nmtypes
  use nmfile
  use nmintegrals
  use nmset
  use nmshell
  use nmharmonic
  use nmbasis
  use nmmatrix
  use nmbessel
  use nmmath
  use nmbecke
  use nmvector

  implicit none

    public :: cfp
    public :: cfp_monomial
    public :: cfp_pack_integrals

    private

      integer :: maxqp
      parameter (maxqp = 1000)

      integer :: llp,nqp,ulp
      integer :: llqp(nmaxshell),ulqp(nmaxshell)
      logical :: skip(nmaxset)
      real(8) :: cfptol,dab,dca,dcb,maxexp,minexp
      real(8) :: rac(3),rbc(3),shlrad(nmaxshell)
      real(8) :: rqp(maxqp),rqp2(maxqp),rqw(maxqp)
      real(8) :: powa(nmaxli*nmaxnco),powb(nmaxli*nmaxnco)

contains

  subroutine cfp(m,pot,basis,h,dbas,cfpfun)
  ! Calculation of Central Force Potential integrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis,pot
    integer :: dbas
    real(8) :: h(dbas,dbas)
    external :: cfpfun
    
    integer :: abas,aset,ashell,bbas,bset,bshell,fa,fb,fsa,fsb
    integer :: ibas,iqp,iset,ishell,jbas,jset,jshell,kset,lmax,na,nb,nqpa
    real(8) :: ra(3),rb(3),rc(3)
    type(nshell) :: sa,sb

    integer :: allocs,dshl,dshp,dpol
    real(8),allocatable :: v(:),rsph(:)
    real(8),allocatable :: shlblk(:,:),ocint(:,:),tblk(:,:)

    ! Initialize
    call cfpini(basis)

    ! Allocate fields 
    lmax = basis_lmax(basis)
    dshp = (2*lmax+1)**2
    dpol = shell_dimgaop(lmax)
    dshl = ((lmax+1)*(lmax+2))/2

    allocate(v(nqp),rsph(dshp),shlblk(dshl,dshl),ocint(dpol,dpol),&
             tblk(dshl,dpol),stat=allocs)
    if (allocs.gt.0) call file_error('CFP: Allocation failure')

    ! Loop over the potential centers 
    do 60 kset=1,pot%nsets
      rc(1:3) = m%atom(pot%set(kset)%atom)%pos(1:3)

      ! Get the potential 
      v(:nqp) = 0.0
      do iqp=1,nqp
        call cfpfun(pot%set(kset),-1,rqp(iqp),v(iqp))
      end do

      ! Screening with respect to this potential 
      call cfpscr(m,pot%set(kset),basis,v,nqpa,'potential')
      if (nqpa.eq.0) go to 60
      call cfpscr(m,pot%set(kset),basis,v,nqpa,'basis')

      ! Weighted potential
      v(1:nqpa) = rqw(1:nqpa)*v(1:nqpa)

      ! Loop on basis functions by atoms 
      fa = 0
      fsa = 0
      do 50 iset=1,basis%nsets
        if (.not.skip(iset)) then
          ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
          ! Geometric parameters for iset
          dca = vector_distance(ra,rc)
          rac(1:3) = ra(1:3) - rc(1:3)
          call cfp_monomial(set_lmax(basis%set(iset)),rac,powa) 
          do 40 ishell=1,basis%set(iset)%nshell
            ashell = fsa + ishell
            sa = basis%set(iset)%shell(ishell)
            na = ((sa%l+1)*(sa%l+2))/2
            if (ulqp(ashell).ge.llqp(ashell)) then
              fb = 0
              fsb = 0
              do 30 jset=1,basis%nsets
                if (.not.skip(jset)) then
                  rb(1:3) = m%atom(basis%set(jset)%atom)%pos(1:3)
                  ! Geometric parameters for jset
                  dab = vector_distance(ra,rb)
                  dcb = vector_distance(rc,rb)
                  rbc(1:3) = rb(1:3) - rc(1:3)
                  call cfp_monomial(set_lmax(basis%set(jset)),rbc,powb) 
                  ! Loop over shells on jset
                  do 20 jshell=1,basis%set(jset)%nshell
                    bshell = fsb + jshell
                    llp = max(llqp(ashell),llqp(bshell))
                    ulp = min(ulqp(ashell),ulqp(bshell))
                    sb = basis%set(jset)%shell(jshell)
                    nb = ((sb%l+1)*(sb%l+2))/2
                    if ((ulp.ge.llp).and.&
                        (shlrad(ashell)+shlrad(bshell).ge.dab).and.&
                        ((jset.gt.iset).or.      &
                         ((jset.eq.iset).and.(jshell.ge.ishell)))) then
                      call cfpblk(0,0,sa,sb,pot%set(kset),ocint,tblk,v,rsph,&
                                  shlblk,dpol,dshl,dshp,cfpfun)
                      do ibas=1,na
                        abas = fa + ibas
                        do jbas=1,nb
                          bbas = fb + jbas
                          h(abas,bbas) = h(abas,bbas) + sa%norm(ibas)*&
                                       sb%norm(jbas)*shlblk(ibas,jbas)
                        end do
                      end do
                    end if
                    fb = fb + nb
 20               continue
                end if
                fsb = fsb + basis%set(jset)%nshell
 30           continue
            end if
            fa = fa + na
 40       continue
        end if
        fsa = fsa + basis%set(iset)%nshell
 50   continue
 60 continue

    deallocate(v,rsph,shlblk,ocint,tblk,stat=allocs)
    if (allocs.gt.0) call file_error('CFP: Deallocation failure')

    ! Symmetrize CFP matrix 
    call matrix_symmetrize(h,dbas,dbas,'uplow')
  end subroutine

  subroutine cfpini(basis)
  ! Central Force Potential integrator INItialization. 
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nbasis) :: basis

    integer :: i,iset,ishell,nshellp
    real(8) :: alpha,d
    type(nshell) :: sa

    cfptol = 1.0e-8
    nqp = 99

    ! Pretabulation of K(z) Bessel functions
    call bessel_initialize

    ! Gauss-Chebyshev abcisas and weights in [0,Infty) interval 
    call becke_radial_quadrature(rqp,rqw,nqp)
    rqp2(1:nqp) = rqp(1:nqp)**2

    minexp = log(ntolnum) - 2.0
    maxexp = 30.0

    ! Initialize screening 
    nshellp = 0
    do iset=1,basis%nsets
      do ishell=1,basis%set(iset)%nshell
        nshellp = nshellp + 1
        sa = basis%set(iset)%shell(ishell)
        ! Find smallest exponent of shell (longest tail)
        alpha = sa%z(1)
        d = abs(sa%c(1))
        do i=2,sa%k
          if (sa%z(i).lt.alpha) then
            alpha = sa%z(i)
            d = abs(sa%c(i))
          end if
        end do
        ! Get the shell radius 
        shlrad(nshellp) = sqrt(shell_gto_radius(alpha,d,sa%l,ntolnum))
      end do
    end do
  end subroutine

  subroutine cfpscr(m,pset,basis,v,nqpa,option)
  ! Central Force Potential integral SCReening.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nmolecule) :: m
    type(nset) :: pset
    type(nbasis) :: basis
    character*(*) :: option
    integer :: nqpa
    real(8) :: v(*)

    integer :: iqp,iset,ishell,nshellp
    real(8) :: rmin,rmax
    real(8) :: ra(3),rc(3)
    type(nshell) :: sa

    ! Potential screening 
    if (option.eq.'potential') then
      nqpa = 0
      do iqp=nqp,1,-1
        if (abs(v(iqp)).gt.ntolnum) then
          nqpa = iqp
          return
        end if
      end do
    ! Basis set screening 
    else if (option.eq.'basis') then
      rc(1:3) = m%atom(pset%atom)%pos(1:3)
      nshellp = 0
      do iset=1,basis%nsets
        ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
        skip(iset) = .true.
        dca = vector_distance(ra,rc)
        do ishell=1,basis%set(iset)%nshell
          nshellp = nshellp + 1
          sa = basis%set(iset)%shell(ishell)

          ! Continuous screening for shells 
          rmin = dca - shlrad(nshellp)
          rmax = dca + shlrad(nshellp)

          iqp = nqpa
          do while ((rqp(iqp).ge.rmin).and.(iqp.gt.1))
            iqp = iqp - 1
          end do
          llqp(nshellp) = iqp + 1

          iqp = nqpa
          do while ((rqp(iqp).gt.rmax).and.(iqp.gt.1))
            iqp = iqp - 1
          end do
          if (rqp(iqp).gt.rmax) then 
            ulqp(nshellp) = 0
          else
            ulqp(nshellp) = iqp 
          end if 
          ! Enable atomic calculation and select first grid 
          ! point (4*i + 1) of the adaptive quadrature 
          if (ulqp(nshellp).ge.llqp(nshellp)) then
            skip(iset) = .false.
            llqp(nshellp) = llqp(nshellp)- mod(llqp(nshellp),4) + 1
          end if 
        end do
      end do
    end if
  end subroutine

  subroutine cfp_monomial(l,r,pow) 
  ! Evaluate geometry factors for Central Force Potential POLynomial expansion.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: l
    real(8) :: r(3),pow(*)
    ! Description of some variables
    ! l  : Maximum angular quantum number of basis set.
    ! pow: Geometric factors from the displacement POWers.
    ! r  : Interatomic distance between basis and potential.

    integer :: lx,ly,lz,ic
    real(8) :: factorx,factory,factorz

    ! Polynomial pre-factors 
    do lx=0,l
      if (lx.eq.0) then
        factorx = 1.0
      else
        factorx = -factorx*r(1)
      end if
      do ly=0,l-lx
        if (ly.eq.0) then
          factory = 1.0
        else
          factory = -factory*r(2)
        end if
        do lz=0,l-lx-ly
          if (lz.eq.0) then
            factorz = 1.0
          else
            factorz = -factorz*r(3)
          end if
          ic = gaop(lx,ly,lz)
          pow(ic) = factorx*factory*factorz
        end do
      end do
    end do
  end subroutine

  subroutine cfpblk(ashift,bshift,sa,sb,pset,ocint,tblk,v,rsph,shlblk,&
                    dpol,dshl,dshp,cfpfun)
  ! Calculation of CFP integral blocks. 
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nshell) :: sa,sb
    type(nset) :: pset
    integer :: ashift,bshift
    integer :: dpol,dshl,dshp
    real(8) :: ocint(dpol,dpol),tblk(dshl,dpol),shlblk(dshl,dshl)
    real(8) :: rsph(dshp),v(*)
    external :: cfpfun

    integer :: i,iqp,j,laa,lab,lbb,ncolaa,ncolbb
    logical :: success
    real(8) :: a2,b2,expo,factor,phi,r,theta,zp
    real(8) :: p(3),sfa(nmaxcon),sfb(nmaxcon),vtmp(maxqp)
    real(8) :: t(0:nmaxli,0:nmaxli)

    ! Initialize 
    laa = sa%l + ashift
    lbb = sb%l + bshift
    lab = laa + lbb
    ocint(1:shell_dimgaop(laa),1:shell_dimgaop(lbb)) = 0.0

    ! Pre-factors 
    a2 = dca**2
    do i=1,sa%k
      sfa(i) = sa%c(i)*exp(-sa%z(i)*a2)
      do j=1,ashift
        sfa(i) = sfa(i)*sa%z(i)
      end do
    end do
    b2 = dcb**2
    do i=1,sb%k
      sfb(i) = sb%c(i)*exp(-sb%z(i)*b2)
      do j=1,bshift
        sfb(i) = sfb(i)*sb%z(i)
      end do
    end do

    ! Loop over primitive Gaussians 
    do i=1,sa%k
      do j=1,sb%k
        ! Angular part 
        zp = -(sa%z(i) + sb%z(j))
        p(:) = 2.0*(sa%z(i)*rac(:) + sb%z(j)*rbc(:))
        call harmonic_xyz_to_spherical(p(1),p(2),p(3),r,theta,phi)
        call harmonic_real_normalized(lab,theta,phi,rsph)
        ! Radial part 
        factor = sfa(i)*sfb(j)
        do iqp=llp,ulp
          expo = zp*rqp2(iqp) + r*rqp(iqp)
          if (expo.gt.maxexp) then
            go to 1000
          else if (expo.lt.minexp) then
            vtmp(iqp) = 0.0
          else
            vtmp(iqp) = factor*v(iqp)*exp(expo)
          end if
        end do
        call cfprad(lab,r,vtmp,t,success)
        if (success) go to 2000 
 1000   continue
        call cfpprim(pset,lab,ashift,bshift,sa%z(i),sb%z(j),sa%c(i),sb%c(j),&
                     r,t,cfpfun)
        ! Build one-center integrals 
 2000   continue
        call cfpocint(laa,lbb,ocint,rsph,t,dpol,dshp)
      end do
    end do
    ! Evaluate the integral using the one-center blocks 
    call cfp_pack_integrals(laa,lbb,ocint,tblk,shlblk,dpol,dshl)
    ! Angular normalization 
    ncolaa = ((laa+1)*(laa+2))/2
    ncolbb = ((lbb+1)*(lbb+2))/2
    do i=1,ncolaa
      do j=1,ncolbb
        shlblk(i,j) = 4.0*pi*shlblk(i,j)
      end do
    end do
  end subroutine

  subroutine cfprad(lmax,wbes,f,t,success)
  ! Central Force Potential RADial quadrature.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: lmax
    logical :: success
    real(8) :: wbes
    real(8) :: t(0:nmaxli,0:nmaxli)
    real(8) :: f(*)

    integer :: qpstep,icyc,iqp,l,llaqp,llqp,n,ulaqp,ulqp
    real(8) :: left,r,right,wr
    real(8) :: besselp(0:nmaxlk),rpow(0:nmaxli)
    real(8) :: tp(0:nmaxli,0:nmaxli),tpp(0:nmaxli,0:nmaxli)

    ! Initialize 
    llqp = 4
    llaqp = llqp + llp - 1 
    ulqp = nqp - llqp + 1
    ulaqp = min(ulp,ulqp)
    qpstep = llqp
    do n=0,lmax
      do l=n,0,-2
        t(l,n) = 0.0 
      end do
    end do
    ! Refinement loop 
    do icyc=1,3 
      ! Quadrature loop 
      do iqp=llaqp,ulaqp,qpstep
        if (abs(f(iqp)).gt.ntolnum) then
          r = rqp(iqp)
          wr = wbes*r
          call bessel_k(lmax,wr,besselp)
          rpow(0) = f(iqp)
          do n=1,lmax
            rpow(n) = rpow(n-1)*r
          end do
          do n=0,lmax
            do l=n,0,-2 
              t(l,n) = t(l,n) + rpow(n)*besselp(l)
            end do
          end do
        end if
      end do
      ! Save current estimates 
      if (icyc.eq.1) then
        do n=0,lmax
          do l=n,0,-2 
            tpp(l,n) = llqp*t(l,n)
          end do
        end do
      else if (icyc.eq.2) then
        do n=0,lmax
          do l=n,0,-2 
            tp(l,n) = llqp*t(l,n)
          end do
        end do
      end if
      ! Next refinement
      if (icyc.lt.3) then
        qpstep = llqp
        llqp = llqp/2
        llaqp = llqp + llp - 1 
        ulqp = nqp - llqp + 1
        ulaqp = min(ulp,ulqp)
      end if
    end do
    ! Check convergence 
    success = .true.
    do n=0,lmax
      do l=n,0,-2 
        left = (t(l,n) - tp(l,n))**2
        right = abs(t(l,n) - tpp(l,n))
        if (left.gt.cfptol*right) then
          if (right.gt.cfptol) then
            success = .false.
            return
          end if
        end if 
      end do
    end do
  end subroutine

  subroutine cfpprim(pset,lmax,ashift,bshift,za,zb,da,db,wbes,t,cfpfun)
  ! Perform one-dimensional adaptive Gauss-Chebyshev quadratures over 
  ! primitive Gaussians for Central Force Potential radial part.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nset) :: pset
    integer :: ashift,bshift,lmax
    real(8) :: da,db,wbes,za,zb
    real(8) :: t(0:nmaxli,0:nmaxli)
    external :: cfpfun

    logical :: converged
    integer :: qpstep,iqp,l,llqp,n,ncycmax,nqpmax,ulqp,nqp_local
    real(8) :: a,b,expf,expo,factor,left,r,right,w
    real(8) :: besselp(0:nmaxlk),rpow(0:nmaxli),wx(maxqp),x(maxqp)
    real(8) :: tp(0:nmaxli,0:nmaxli),tpp(0:nmaxli,0:nmaxli),tw(0:nmaxli,0:nmaxli)

    ! Calculate linearly mapped quadrature points 
    nqpmax = nqp
    ncycmax = 0
    do while (nqpmax.lt.maxqp)
      nqpmax = 2*nqpmax + 1
      ncycmax = ncycmax + 1
    end do
    nqpmax = (nqpmax-1)/2
    ncycmax = ncycmax + 2
    call becke_tskgc(x,wx,nqpmax)
    ! Linear mapping 
    call cfpmap(za,zb,dca,dcb,a,b,'initialize')
    ! Initialize quadrature 
    factor = da*db
    if (ashift.gt.0) factor = factor*za**ashift
    if (bshift.gt.0) factor = factor*zb**bshift
    converged = .false.
    nqp_local = ((nqp - 1)/2 - 1)/2
    llqp = 2**(ncycmax - 1)
    ulqp = nqpmax - llqp + 1
    qpstep = llqp
    do n=0,lmax
      do l=n,0,-2 
        tw(l,n) = 0.0 
      end do
    end do
    ! Quadrature loop 
    do while (.not.converged) 
      do iqp=llqp,ulqp,qpstep
        ! Calculate abscissa and weight 
        call cfpmap(a,b,x(iqp),wx(iqp),r,w,'mapping')
        expo = - za*(dca-r)**2 - zb*(dcb-r)**2
        if (expo.ge.minexp) then
          ! Calculate exponential function on the grid 
          call cfpfun(pset,-1,r,expf)
          expf = w*factor*expf*exp(expo)
          ! Calculate radial powers on the grid 
          rpow(0) = expf
          do l=1,lmax
            rpow(l) = rpow(l-1)*r
          end do
          ! Calculate Bessel functions
          call bessel_k(lmax,wbes*r,besselp)
          ! Integral estimates
          do n=0,lmax
            do l=n,0,-2 
              tw(l,n) = tw(l,n) + rpow(n)*besselp(l)
            end do
          end do
        end if
      end do
      ! Save integral estimates 
      do n=0,lmax
        do l=n,0,-2 
          tpp(l,n) = tp(l,n)
          tp(l,n) = t(l,n)
          t(l,n) = llqp*tw(l,n)
        end do
      end do
      ! Test convergence 
      if (nqp_local.ge.nqp) then
        converged = .true.
        do n=0,lmax
          do l=n,0,-2 
            left = (t(l,n) - tp(l,n))**2
            right = abs(t(l,n) - tpp(l,n))
            if (left.gt.cfptol*right) then
              if (right.gt.cfptol) then
                converged = .false.
                go to 1000
              end if
            end if 
          end do
        end do
 1000   continue
        if ((nqp_local.eq.nqpmax).and.(.not.converged)) then
          call file_error('CFPPRIM: Allocation failure')
        end if
      end if
      ! Prepare next cycle 
      nqp_local = 2*nqp_local + 1
      qpstep = llqp
      llqp = llqp/2
      ulqp = nqpmax - llqp + 1
    end do
  end subroutine


  subroutine cfpocint(la,lb,ocint,rsph,t,dpol,dshp)
  ! Calculation of CFP One-Center INTegral blocks. 
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: dpol,dshp
    real(8) :: ocint(dpol,dpol),rsph(dshp),t(0:nmaxli,0:nmaxli)

    integer :: kax,kay,kaz,kbx,kby,kbz,kx,ky,kz
    integer :: fl,k,kac,kbc,kc,l,la,lb,ish
    real(8) :: factor

    ! Left loop 
    do kax=0,la
      do kay=0,la-kax
        do kaz=0,la-kax-kay
          kac = gaop(kax,kay,kaz)
          ! Right loop 
          do kbx=0,lb
            kx = kax + kbx
            do kby=0,lb-kbx
              ky = kay + kby
              do kbz=0,lb-kbx-kby 
                kz = kaz + kbz
                kbc = gaop(kbx,kby,kbz)
                ! Build one-center integrals 
                kc = gaop(kx,ky,kz)
                k = kx + ky + kz 
                do l=k,0,-2 
                  fl = l*l
                  factor = rsph(fl+1)*ptostm(1,kc,l)
                  do ish=2,2*l+1
                    factor = factor + rsph(fl+ish)*ptostm(ish,kc,l)
                  end do
                  ocint(kac,kbc) = ocint(kac,kbc) + factor*t(l,k)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine

  subroutine cfp_pack_integrals(la,lb,ocint,tblk,shlblk,dpol,dshl)
  ! Evaluate CFP INTegrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: dpol,dshl,la,lb
    real(8) :: ocint(dpol,dpol),tblk(dshl,dpol),shlblk(dshl,dshl)

    integer :: ax,ay,az,kax,kay,kaz,as,kac,akc,ckax,ckay,ckaz
    integer :: bx,by,bz,kbx,kby,kbz,bs,kbc,bkc,ckbx,ckby,ckbz
    integer :: ncola,ncosetlb
    real(8) :: bax,bay,baz
    real(8) :: bbx,bby,bbz

    ! Initialize 
    shlblk(:,:) = 0.0
    ncola = ((la+1)*(la+2))/2
    ncosetlb = shell_dimgaop(lb)
    tblk(1:ncola,1:ncosetlb) = 0.0

    ! Left side expansion 
    do ax=0,la
      do ay=0,la-ax
        az = la - ax - ay
        as = raop(ax,ay,az)
        do kax=0,ax
          bax = math_binomial(ax,kax)
          ckax = ax - kax
          do kay=0,ay
            bay = bax*math_binomial(ay,kay)
            ckay = ay - kay
            do kaz=0,az
              ckaz = az - kaz
              akc = gaop(ckax,ckay,ckaz)
              baz = bay*powa(akc)*math_binomial(az,kaz)
              if (abs(baz).gt.ntolnum) then
                kac = gaop(kax,kay,kaz)
                do kbc=1,ncosetlb
                  tblk(as,kbc) = tblk(as,kbc) + baz*ocint(kac,kbc)
                end do
              end if
            end do
          end do
        end do
      end do
    end do

    ! Right side expansion 
    do bx=0,lb
      do by=0,lb-bx
        bz = lb - bx - by
        bs = raop(bx,by,bz)
        do kbx=0,bx
          bbx = math_binomial(bx,kbx)
          ckbx = bx - kbx
          do kby=0,by
            bby = bbx*math_binomial(by,kby)
            ckby = by - kby
            do kbz=0,bz
              ckbz = bz - kbz
              bkc = gaop(ckbx,ckby,ckbz)
              bbz = bby*math_binomial(bz,kbz)*powb(bkc)
              if (abs(bbz).gt.ntolnum) then
                kbc = gaop(kbx,kby,kbz)
                do as=1,ncola
                  shlblk(as,bs) = shlblk(as,bs) + tblk(as,kbc)*bbz
                end do
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine

  subroutine cfpmap(inp1,inp2,inp3,inp4,out1,out2,option)
  ! Central Force Potential quadrature MAPping utility.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    !
    ! option: a) initialize: Initalize CFP mapping.
    !         b) mapping: Map quadrature on primitives.
    !
    ! For option = initialize: 
    !   inp1: First basis set exponent.
    !   inp2: Second basis set exponent.
    !   inp3: Interatomic distance from potential to first basis centers.
    !   inp4: Interatomic distance from potential to second basis centers.
    !   out1: First mapping parameter (slope).
    !   out1: Second mapping parameter (origin ordinate).
    !
    ! For option = mapping: 
    !   inp1: First mapping parameter (slope).
    !   inp1: Second mapping parameter (origin ordinate).
    !   inp3: Unmapped abcissa. 
    !   inp4: Unmapped weight. 
    !   out1: Mapped abcissa. 
    !   out1: Mapped weight. 

    character*(*) :: option
    real(8) :: inp1,inp2,inp3,inp4,out1,out2

    real(8) :: p,rmax,rmin,sigma,zp

    ! Linear mapping: (-1,1) ==> (r_min,r_max) 
    if (option.eq.'initialize') then
      ! Get parameters 
      zp = inp1 + inp2
      p = (inp1*inp3 + inp2*inp4)/zp
      sigma = 1.0/sqrt(zp)
      ! Radial limits 
      rmin = max(0.0,p - 7.0*sigma)
      rmax = p + 9.0*sigma
      ! Mapping factors 
      out1 = 0.5*(rmax - rmin)
      out2 = 0.5*(rmax + rmin)
      ! Do the mapping 
    else if (option.eq.'mapping') then
      out1 = inp2 + inp1*inp3 
      out2 = inp1*inp4
    end if
  end subroutine

end module
