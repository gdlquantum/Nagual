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
module nmecp
! Literature:
! Half‐numerical evaluation of pseudopotential integrals
! R. Flores‐Moreno et al. Journal of computational chemistry 27 (9), 1009-1019

  use nmtypes
  use nmbasis
  use nmfile
  use nmintegrals
  use nmmath
  use nmset
  use nmshell
  use nmmatrix
  use nmbessel
  use nmharmonic
  use nmcfp
  use nmbecke
  use nmvector

  implicit none

    public :: ecp

    private
      integer :: maxqp
      parameter (maxqp = 1000)

      integer :: llp,nqp,ulp
      integer :: llqp(nmaxshell),ulqp(nmaxshell)
      logical :: skip(nmaxset)
      real(8) :: ecptol,dca,dcb,maxexp,minexp
      real(8) :: rac(3),rbc(3),shlrad(nmaxshell)
      real(8) :: rqp(maxqp),rqp2(maxqp),rqw(maxqp)
      real(8) :: powa(nmaxli*nmaxnco),powb(nmaxli*nmaxnco)
      real(8) :: rint(0:nmaxli,0:nmaxli,0:nmaxli)         

      ! Description of module variables
      !
      ! ecptab  : Table of operator part for radial quadratures.
      ! h       : Core Hamiltonian matrix.

contains

  subroutine ecp(m,pot,basis,h,dbas)
  ! Calculate ECP integrals and add them to the core Hamiltonian matrix.
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis,pot
    integer :: dbas
    real(8) :: h(dbas,dbas)

    ! Calculate local ECP integrals 
    call cfp(m,pot,basis,h,dbas,ecppotf)

    ! Calculate semi-local ECP integrals
    call semecp(m,pot,basis,h,dbas)
  end subroutine

  subroutine ecppotf(pset,lc,r,pot)
  ! Calculate ECP Central Force Potential.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nset) :: pset
    integer :: lc
    real(8) :: pot,r

    integer :: i,ishell,l,n
    real(8) :: r2,radial,expf
    type(nshell) :: sa

    pot = 0.0
    ! Quick return, if no ECP is requested 
    if (pset%nshell.le.0) return
    ! If local, set L
    if (lc.lt.0) then
      l = set_lmax(pset) 
    else
      l = lc
    end if
    ! Pre-compute r^2 
    r2 = r*r
    ! Loop over all ECP shells of the current atom
    do ishell=1,pset%nshell
      sa = pset%shell(ishell)
      if (sa%l.eq.l) then  
        ! Calcualte the radial prefactor of the potential
        n = sa%n + 2
        if (n.eq.0) then
          radial = 1.0
        else 
          radial = r**n
        end if
        ! Contract the potential of the current shell
        expf = 0.0
        do i=1,sa%k
          expf = expf + sa%c(i)*exp(-sa%z(i)*r2)
        end do
        pot = pot + radial*expf
      end if
    end do
  end subroutine

  subroutine semecp(m,pot,basis,h,dbas)
  ! Calculation of semi-local ECP integrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis,pot
    integer :: dbas
    real(8) :: h(dbas,dbas)

    integer :: abas,ashell,bbas,bshell,fa,fb,fsa,fsb
    integer :: ibas,iset,ishell,jbas,jset,jshell,kset
    integer :: lac,lamax,lbc,lbmax,lcmax,na,nb,nqpa,nshl
    real(8) :: ra(3),rb(3),rc(3)
    type(nshell) :: sa,sb

    integer :: allocs,decp,diml,dmu,dpoc,dpol,dshl,dshp
    real(8),allocatable :: angint(:,:,:),ecptab(:,:,:),bastab(:,:,:),&
                           ocint(:,:),shlblk(:,:),taba(:,:),tabb(:,:),&
                           tblk(:,:),slaa(:,:,:),slab(:,:,:)

    call ecpini(basis,nshl)

    ! dimensions 
    decp = basis_lmax(pot) - 1
    dmu = basis_lmax(basis)
    diml = dmu + max(decp,dmu) 
    dpoc = shell_dimgaop(diml)
    dpol = shell_dimgaop(dmu)
    dshl = ((dmu+1)*(dmu+2))/2
    dshp = harmonic_pointer_size(dmu+decp)

    allocate(angint(dshp,dshp,dpoc),ecptab(nqp,0:decp,0:diml),&
      bastab(nqp,nshl,0:diml),ocint(dpol,dpol),tblk(dshl,dpol),&
      shlblk(dshl,dshl),taba(nqp,0:diml),tabb(nqp,0:diml),&
      slaa(0:diml,dshp,dpol),slab(0:diml,dshp,dpol),stat=allocs)
    if (allocs.gt.0) call file_error('ECP: Allocation failure')

    ! Angular integrals 
    call ecpang(0,dmu,decp,angint,dshp,dpoc)

    do 50 kset=1,pot%nsets
      if (pot%set(kset)%nshell.eq.0) go to 50
      rc(1:3) = m%atom(pot%set(kset)%atom)%pos(1:3)
      lcmax = set_lmax(pot%set(kset)) - 1 
      ! Pre-tabulation 
      call ecppot(m,pot%set(kset),ecptab,nqp,decp,diml,nqpa)
      if (nqpa.eq.0) go to 50
      call ecpscr(m,pot%set(kset),basis,ecptab(:,0,0),nqpa,'basis')
      call ecpbas(0,.true.,m,pot%set(kset),basis,bastab,nshl,diml)
      fa = 0
      fsa = 0
      do 40 iset=1,basis%nsets
        if (.not.skip(iset)) then
          ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
          lamax = set_lmax(basis%set(iset)) 
          dca = vector_distance(ra,rc)
          rac(1:3) = ra(1:3) - rc(1:3)
          call ecpgeo(lamax,lcmax,rac,dca,angint,powa,slaa,diml,dshp,dpoc,dpol) 
          do 30 ishell=1,basis%set(iset)%nshell
            ashell = fsa + ishell
            sa = basis%set(iset)%shell(ishell)
            na = ((sa%l+1)*(sa%l+2))/2
            llp = llqp(ashell)
            ulp = ulqp(ashell)
            if (ulp.ge.llp) then
              lac = sa%l + lcmax 
              taba(llp:ulp,0:lac) = bastab(llp:ulp,ashell,0:lac)
              fb = 0
              fsb = 0
              do 20 jset=1,basis%nsets
                if (.not.skip(jset)) then
                  rb(1:3) = m%atom(basis%set(jset)%atom)%pos(1:3)
                  lbmax = set_lmax(basis%set(jset))
                  dcb = vector_distance(rb,rc)
                  rbc(1:3) = rb(1:3) - rc(1:3)
                  call ecpgeo(lbmax,lcmax,rbc,dcb,angint,powb,slab,diml,dshp,dpoc,dpol)
                  do 10 jshell=1,basis%set(jset)%nshell
                    bshell = fsb + jshell
                    sb = basis%set(jset)%shell(jshell)
                    nb = ((sb%l+1)*(sb%l+2))/2
                    llp = max(llqp(ashell),llqp(bshell))
                    ulp = min(ulqp(ashell),ulqp(bshell))
                    if ((ulp.ge.llp).and.((jset.gt.iset).or.&
                         ((jset.eq.iset).and.(jshell.ge.ishell)))) then
                      lbc = sb%l + lcmax
                      tabb(llp:ulp,0:lbc) = bastab(llp:ulp,bshell,0:lbc)
                      call ecpblk(0,0,sa,sb,pot%set(kset),taba,tabb,ecptab,&
                                  slaa,slab,ocint,tblk,shlblk,diml,decp,dshp,&
                                  dpol,dshl)
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
 10               continue
                end if
                fsb = fsb + basis%set(jset)%nshell
 20           continue
            end if
            fa = fa + na
 30       continue
        end if
        fsa = fsa + basis%set(iset)%nshell
 40   continue
 50 continue

    deallocate(angint,ecptab,bastab,ocint,tblk,shlblk,taba,tabb,&
      slaa,slab,stat=allocs)
    if (allocs.gt.0) call file_error('ECP: Deallocation failure')

    ! Symmetrize ECP matrix 
    call matrix_symmetrize(h,dbas,dbas,'uplow')
  end subroutine 

  subroutine ecpang(l,lbasmax,lecpmax,angint,dshp,dpoc)
  ! ECP semi-local angular integrals.
  ! < S(LA,MA) | S(LB,MB) | x^i y^j z^k >  
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: dpoc,dshp,l,lbasmax,lecpmax
    real(8) :: angint(dshp,dshp,dpoc)
    ! Description of some variables
    ! dpoc: Dimension of the coupled polynomials.
    ! dshp: Dimension of the spherical harmonics. 
    ! angint: Analytic angular integrals for ECPs.

    integer :: la,lb,lmax,lp,ma,mb,nsh,par,px,py,pz,sha,shb

    ! Pre-tabulation of the angular integrals 
    lmax = lbasmax + lecpmax + l 
    nsh = harmonic_pointer_size(lmax)
    angint(1:nsh,1:nsh,:) = 0.0

    do la=0,lmax
      do lb=0,lmax
        par = max(mod(la+lb,2),la-lb)
        do lp=par,lbasmax+l,2
          do ma=1,2*la+1
            sha = la*la+ma 
            do mb=1,2*lb+1
              shb = lb*lb+mb 
              do px=0,lp
                do py=0,lp-px
                  pz = lp - px - py
                  angint(sha,shb,gaop(px,py,pz)) = &
                    harmonic_ssxyz_integral(la,ma,lb,mb,px,py,pz)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine

  subroutine ecppot(m,pset,tab,nqp,decp,diml,nqpa)
  ! Pre-tabulation of effective core potential functions for 
  ! the numerical integration of the radial ECP integrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nmolecule) :: m
    type(nset) :: pset
    integer :: decp,diml,nqp,nqpa
    real(8) :: tab(nqp,0:decp,0:diml)

    integer :: lcmax,llqp,ulqp

    integer :: iqp,lc,n,nqpashl
    real(8) :: v(maxqp)
    type(nbasis) :: dummy

    ! Initialization 
    nqpa = 0
    lcmax = set_lmax(pset) - 1
    tab(1:nqp,0:lcmax,0:diml) = 0.0 

    ! Set grid point limits 
    llqp = 1
    ulqp = nqp
    ! Loop over all (semi-local) ECP shells 
    do 10 lc=0,lcmax
      v(1:nqp) = 0.0
      do iqp=llqp,ulqp
        call ecppotf(pset,lc,rqp(iqp),v(iqp))
      end do
      ! Potential screening 
      call ecpscr(m,pset,dummy,v,nqpashl,'potential')
      nqpa = max(nqpashl,nqpa)
      ! Weighting 
      tab(nqpashl+1:nqp,lc,0) = 0.0
      if (nqpashl.eq.0) go to 10
      tab(1:nqpashl,lc,0) = rqw(1:nqpashl)*v(1:nqpashl)
      ! Generate higher L entries 
      do n=1,diml
        tab(1:nqpashl,lc,n) = rqp(1:nqpashl)*tab(1:nqpashl,lc,n-1)
      end do
 10 continue
  end subroutine

  subroutine ecpbas(n,scaled,m,pset,basis,tab,dshl,diml)
  ! Pre-tabulation of basis functions for the numerical
  ! integration of the radial ECP integrals.
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: basis
    type(nset) :: pset
    integer :: dshl,diml,n
    logical :: scaled
    real(8) :: tab(nqp,dshl,0:diml)

    integer :: ashell,i,iqp,iset,ishell,j,la,lmax
    real(8) :: expo,factor,wz,z
    real(8) :: ra(3),rc(3),bessel(0:nmaxlk)
    type(nshell) :: sa

    ! Pre-tabulation of basis functions 
    tab(1:nqp,1:dshl,0:diml) = 0.0
 
    rc(1:3) = m%atom(pset%atom)%pos(1:3)

    ! Loop over all atoms 
    ashell = 0
    do 50 iset=1,basis%nsets
      if (skip(iset)) go to 50
      ra(1:3) = m%atom(basis%set(iset)%atom)%pos(1:3)
      dca = vector_distance(ra,rc)
      ! Loop over all shells of the current atom 
      do ishell=1,basis%set(iset)%nshell
        ashell = ashell + 1
        sa = basis%set(iset)%shell(ishell)
        lmax = sa%l + n + set_lmax(pset) - 1
        ! Loop over all primitives of the current shell 
        do i=1,sa%k
          z = sa%z(i)
          factor = sa%c(i)
          do j=1,n
            if (scaled) factor = factor*z
          end do
          ! Do the grid work 
          do iqp=llqp(ashell),ulqp(ashell)
            wz = 2.0*z*dca*rqp(iqp)
            call bessel_k(lmax,wz,bessel)
            expo = exp(-z*(dca - rqp(iqp))**2)
            bessel(0:lmax) = bessel(0:lmax)*expo
            tab(iqp,ashell,0:lmax) = tab(iqp,ashell,0:lmax) + &
            factor*bessel(0:lmax)
          end do
        end do
      end do
 50 continue
  end subroutine

  subroutine ecpgeo(la,lc,a,r,angint,pow,slap,diml,dshp,dpoc,dpol) 
  ! Evaluate factors which depend only on the geometry.
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    integer :: diml,dshp,dpoc,dpol
    real(8) :: angint(dshp,dshp,dpoc),slap(0:diml,dshp,dpol)
    real(8) :: a(3),pow(*)
    integer :: isph,fl,l,la,lc,lmax,m,npol,nsphc
    real(8) :: phi,r,theta
    real(8) :: apc(dpoc),rsph((nmaxli+1)**2)
    integer :: dimcop,dimshp  

    ! Loop limits
    lmax = lc + la 
    nsphc = (lc+1)**2
    npol = shell_dimgaop(la)

    ! Spherical coordinates 
    call harmonic_xyz_to_spherical(a(1),a(2),a(3),r,theta,phi)

    ! Polynomial prefactors 
    call ecppol(la,a,pow) 

    ! Evaluate normalized real spherical harmonics 
    call harmonic_real_normalized(lmax,theta,phi,rsph)

    ! Build angular blocks 
    do isph=1,nsphc
      do l=0,lmax
        fl = l*l
        apc(1:npol) = 0.0
        do m=1,2*l+1
          apc(1:npol) = apc(1:npol) + rsph(fl+m)*angint(fl+m,isph,1:npol)     
        end do
        slap(l,isph,1:npol) = apc(1:npol)
      end do
    end do
  end subroutine


  subroutine ecpblk(ashift,bshift,sa,sb,pset,a,b,c,slaa,slab,ocint,tblk,&
                        shlblk,diml,decp,dshp,dpol,dshl)
  ! Add the contribution of the ECP semi-local integrals.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nset) :: pset
    type(nshell) :: sa,sb
    integer :: ashift,bshift

    integer :: diml,decp,dshp,dpol,dshl
    real(8) :: a(nqp,0:diml),b(nqp,0:diml),c(nqp,0:decp,0:diml)
    real(8) :: slaa(0:diml,dshp,dpol),slab(0:diml,dshp,dpol)
    real(8) :: ocint(dpol,dpol),tblk(dshl,dpol),shlblk(dshl,dshl)
    logical :: success
    integer :: i,j,laa,lbb,lc,patom
    integer :: ncolaa,ncolbb,ncoslaa,ncoslbb

    ! Initialize 
    laa = sa%l + ashift
    lbb = sb%l + bshift
    ncoslaa = shell_dimgaop(laa)
    ncoslbb = shell_dimgaop(lbb)
    ocint(1:ncoslaa,1:ncoslbb) = 0.0
    ! Compute ECP one-center integrals
    do lc=0,set_lmax(pset)-1
      ! Use pre-tabulation over contracted Gaussians 
      call ecprad(laa,lbb,lc,a,b,c,nqp,diml,decp,success)
      ! If necessary, integrate over primitive Gaussians 
      if (.not.success) then
        call ecpprim(ashift,bshift,pset,sa,sb,lc)
      end if
      ! Save the one-center integral block 
      call ecpocint(laa,lbb,lc,slaa,slab,ocint,diml,dpol,dshp)
    end do
    ! Evaluate the integral using the one-center blocks 
    call ecpint(laa,lbb,ocint,tblk,shlblk,dpol,dshl)
    ! Angular normalization 16 pi^2 
    ncolaa = ((laa+1)*(laa+2))/2
    ncolbb = ((lbb+1)*(lbb+2))/2
    do i=1,ncolaa
      do j=1,ncolbb
        shlblk(i,j) = 16.0*pi**2*shlblk(i,j)
      end do
    end do
  end subroutine

  subroutine ecprad(la,lb,lc,a,b,c,nqp,diml,decp,success)
  ! Perform adaptive Gauss-Chebyshev quadratures for the 
  ! semi-local radial part using pre-tabulated functions.
  ! Roberto Flores-Moreno, 2004, 2018
  implicit none
    integer :: decp,diml,nqp
    real(8) :: a(nqp,0:diml),b(nqp,0:diml),c(nqp,0:decp,0:diml)
    logical :: success
    integer :: la,laa,lab,lac,lb,lbb,lbc,lc,l,iqp,even,odd
    real(8) :: qt(0:nmaxl),qtp(0:nmaxl),qtpp(0:nmaxl)
    real(8) :: oddpot(0:nmaxl,maxqp),evenpot(0:nmaxl,maxqp)

    ! Limits and initialization 
    lab = la + lb
    lac = la + lc
    lbc = lb + lc
    even = lab/2
    odd = (lab-1)/2
    success = .true.
    do l=0,even
      if (l+l+1.le.diml) then
        oddpot(l,llp:ulp) = c(llp:ulp,lc,l+l+1)
      end if
      evenpot(l,llp:ulp) = c(llp:ulp,lc,l+l)
    end do
    ! Gauss-Chebyshev quadrature with pre-tabulated functions 
    do laa=0,lac
      do lbb=0,lbc
        ! Radial integrals with even radial power 
        if (mod(laa+lbb,2).eq.0) then
          ! Initialize quadrature table
          qt(0:even) = 0.0
          ! Get the (p-1)/2 points quadrature estimates 
          do iqp=llp+3,ulp,4
            qt(0:even) = qt(0:even) + a(iqp,laa)*b(iqp,lbb)* &
                         evenpot(0:even,iqp)
          end do
          qtpp(0:even) = 4.0*qt(0:even)
          ! Get the p points quadrature estimates 
          do iqp=llp+1,ulp,4
            qt(0:even) = qt(0:even) + a(iqp,laa)*b(iqp,lbb)* &
                         evenpot(0:even,iqp)
          end do
          qtp(0:even) = 2.0*qt(0:even)
          ! Get the 2p+1 points quadrature estimates 
          do iqp=llp,ulp,2
            qt(0:even) = qt(0:even) + a(iqp,laa)*b(iqp,lbb)* &
                         evenpot(0:even,iqp)
          end do
          ! Check convergence and save integrals 
          qtp(0:even) = (qt(0:even) - qtp(0:even))**2 - &
                        ecptol*abs(qt(0:even) - qtpp(0:even))
          do l=0,even
            if (qtp(l).gt.0.0) then
              if (abs(qt(l)-qtpp(l)).gt.ecptol) then
                success = .false.
                return 
              end if  
            end if
            rint(l+l,laa,lbb) = qt(l)
          end do
        ! Radial integrals with odd radial power 
        else if (odd.ge.0) then
          ! Initialize quadrature table 
          qt(0:odd) = 0.0
          ! Get the (p-1)/2 points quadrature estimates 
          do iqp=llp+3,ulp,4
            qt(0:odd) = qt(0:odd) + a(iqp,laa)*b(iqp,lbb)*&
                        oddpot(0:odd,iqp)
          end do
          qtpp(0:odd) = 4.0*qt(0:odd)
          ! Get the p points quadrature estimates 
          do iqp=llp+1,ulp,4
            qt(0:odd) = qt(0:odd) + a(iqp,laa)*b(iqp,lbb)*&
                        oddpot(0:odd,iqp)
          end do
          qtp(0:odd) = 2.0*qt(0:odd)
          ! Get the 2p+1 points quadrature estimates
          do iqp=llp,ulp,2
            qt(0:odd) = qt(0:odd) + a(iqp,laa)*b(iqp,lbb)* &
                        oddpot(0:odd,iqp)
          end do
          ! Check convergence and save integrals 
          qtp(0:odd) = (qt(0:odd) - qtp(0:odd))**2 - &
                       ecptol*abs(qt(0:odd) - qtpp(0:odd))
          do l=0,odd
            if (qtp(l).gt.0.0) then
              if (abs(qt(l)-qtpp(l)).gt.ecptol) then  
                success = .false.
                return 
              end if  
            end if
            rint(l+l+1,laa,lbb) = qt(l)
          end do
        end if
      end do
    end do
  end subroutine

  subroutine ecpprim(ashift,bshift,pset,sa,sb,lc)
  ! Perform adaptive Gauss-Chebyshev quadratures for the semilocal radial part.
  ! Roberto Flores-Moreno, 2003, 2018
  implicit none
    integer :: ashift,bshift
    type(nset) :: pset
    type(nshell) :: sa,sb

    logical :: converged
    integer :: qpstep
    integer :: laa,lab,lac,lbb,lbc,lc
    integer :: i,iqp,j,l,llqp,nqp_local,nqpmax,ncycmax,ulqp
    real(8) :: a,b,expf,expo,factor,ka,kb,left,r,right,w,za,zb
    real(8) :: bessela(0:nmaxlk),besselb(0:nmaxlk),rpow(0:nmaxli),&
         wx(maxqp),x(maxqp)

    real(8) :: qt(0:2*nmaxl,0:nmaxli,0:nmaxli)
    real(8) :: qtp(0:2*nmaxl,0:nmaxli,0:nmaxli)
    real(8) :: qtpp(0:2*nmaxl,0:nmaxli,0:nmaxli)
    real(8) :: qtw(0:2*nmaxl,0:nmaxli,0:nmaxli)

    ! Limits and initialization 
    lab = sa%l + sb%l
    lac = sa%l + lc
    lbc = sb%l + lc
    rint(0:lab,0:lac,0:lbc) = 0.0
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
    ! Primitive loop 
    do i=1,sa%k
      do j=1,sb%k
        ! Linear mapping 
        za = sa%z(i)
        zb = sb%z(j)
        call ecpmap(za,zb,dca,dcb,a,b,'initialize')
        ! Weights 
        ka = 2.0*dca*za
        kb = 2.0*dcb*zb
        factor = sa%c(i)*sb%c(j)
        do l=1,ashift
          factor = factor*za
        end do
        do l=1,bshift
          factor = factor*zb
        end do
        ! Initialize quadrature 
        converged = .false.
        nqp_local = ((nqp - 1)/2 - 1)/2
        llqp = 2**(ncycmax-1)
        ulqp = nqpmax - llqp + 1
        qpstep = llqp
        qtw(0:lab,0:lac,0:lbc) = 0.0
        ! Quadrature loop 
        do while (.not.converged) 
          do iqp=llqp,ulqp,qpstep
            ! Get mapped abcisa and weight 
            call ecpmap(a,b,x(iqp),wx(iqp),r,w,'mapping')
            expo = - za*(dca-r)**2 - zb*(dcb-r)**2
            if (expo.ge.minexp) then
              ! Calculate Bessel functions 
              call bessel_k(lac,ka*r,bessela)
              call bessel_k(lbc,kb*r,besselb)
              ! Calculate potential function on the grid 
              call ecppotf(pset,lc,r,expf) 
              expf = w*exp(expo)*expf
              ! Calculate radial powers on the grid 
              rpow(0) = factor*expf
              do l=1,lab
                rpow(l) = rpow(l-1)*r
              end do
              ! Integral estimates 
              do laa=0,lac
                do lbb=0,lbc
                  qtw(0:lab,laa,lbb) = qtw(0:lab,laa,lbb) + &
                  rpow(0:lab)*bessela(laa)*besselb(lbb)
                end do
              end do
            end if
          end do
          qtpp(0:lab,0:lac,0:lbc) = qtp(0:lab,0:lac,0:lbc)
          qtp(0:lab,0:lac,0:lbc) = qt(0:lab,0:lac,0:lbc)
          qt(0:lab,0:lac,0:lbc) = llqp*qtw(0:lab,0:lac,0:lbc)
          ! Convergence test
          if (nqp_local.ge.nqp) then
            converged = .true.
            do l=0,lab
              do laa=0,lac
                do lbb=0,lbc
                  left = (qt(l,laa,lbb) - qtp(l,laa,lbb))**2
                  right = abs(qt(l,laa,lbb) - qtpp(l,laa,lbb))
                  if (left.gt.ecptol*right) then
                    if (right.gt.ecptol) then
                      converged = .false.
                      go to 1000
                    end if
                  end if 
                end do
              end do
            end do
 1000       continue
            if ((nqp_local.eq.nqpmax).and.(.not.converged)) then
              call file_error('ecpprim: ecp quadrature failed')
            end if
          end if
          ! Prepare next cycle 
          nqp_local  = 2*nqp_local + 1
          qpstep = llqp
          llqp = llqp/2
          ulqp = nqpmax - llqp + 1
        end do
        ! Store integral values
        rint(0:lab,0:lac,0:lbc) = rint(0:lab,0:lac,0:lbc)+qt(0:lab,0:lac,0:lbc)
      end do
    end do
  end subroutine

  subroutine ecpocint(la,lb,lc,slaa,slab,ocint,diml,dpol,dshp)
  ! Evaluate ECP one-center semi-local integral terms.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    integer :: diml,dpol,dshp
    real(8) :: slaa(0:diml,dshp,dpol),slab(0:diml,dshp,dpol),ocint(dpol,dpol)
    integer :: la,laa,lac,lacmax,lab,lb,lbb,lbc,lbcmax,lc,llac,llbc
    integer :: l1,l2,llsh,para,parb,ulsh
    real(8) :: termint

    ! Limits 
    lacmax = la + lc
    lbcmax = lb + lc
    llsh = lc*lc + 1
    ulsh = llsh + 2*lc
    ! Loop over polynoms derived from center A 
    do laa=0,la
      if (laa.eq.0) then
        llac = 1
        para = lc
      else 
        llac = shell_dimgaop(laa-1) + 1
        para = max(mod(laa+lc,2),lc-laa)
      end if
      do lac=llac,shell_dimgaop(laa)
        ! Loop over polynoms derived from center B 
        do lbb=0,lb
          lab = laa + lbb
          if (lbb.eq.0) then
            llbc = 1
            parb = lc
          else 
            llbc = shell_dimgaop(lbb-1) + 1
            parb = max(mod(lc+lbb,2),lc-lbb)
          end if
          do lbc=llbc,shell_dimgaop(lbb)
            ! Initialize integral value 
            termint = 0.0
            ! Evaluate the integral 
            do l1=para,lacmax,2
              do l2=parb,lbcmax,2
                termint = termint + rint(lab,l1,l2)*&
                sum(slaa(l1,llsh:ulsh,lac)*slab(l2,llsh:ulsh,lbc))
              end do
            end do
            ! Update one-center integral array 
            ocint(lac,lbc) = ocint(lac,lbc) + termint
          end do
        end do
      end do
    end do
  end subroutine

  subroutine ecpscr(m,pset,b,v,nqpa,option)
  ! ECP integral SCReening.
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nmolecule) :: m
    type(nset) :: pset
    type(nbasis) :: b
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
      do iset=1,b%nsets
        ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
        skip(iset) = .true.
        dca = vector_distance(ra,rc)
        do ishell=1,b%set(iset)%nshell
          nshellp = nshellp + 1
          sa = b%set(iset)%shell(ishell)

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

  subroutine ecpmap(inp1,inp2,inp3,inp4,out1,out2,option)
  ! ECP quadrature MAPping utility.
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
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

  subroutine ecppol(l,r,pow) 
  ! Evaluate geometry factors for ECP POLynomial expansion.
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

  subroutine ecpint(la,lb,ocint,tblk,shlblk,dpol,dshl)
  ! Evaluate ECP INTegrals.
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

  subroutine ecpini(basis,dshl)
  ! ECP integrator INItialization. 
  ! Roberto Flores-Moreno, 2005, 2018
  implicit none
    type(nbasis) :: basis
    integer :: dshl

    integer :: i,iset,ishell
    real(8) :: alpha,d
    type(nshell) :: sa

    ! Pretabulation of K(z) Bessel functions
    call bessel_initialize

    ecptol = 1.0e-8
    nqp = 99

    ! Gauss-Chebyshev abcisas and weights in [0,Infty) interval 
    call becke_radial_quadrature(rqp,rqw,nqp)
    rqp2(1:nqp) = rqp(1:nqp)**2

    minexp = log(ntolnum) - 2.0
    maxexp = 30.0

    ! Initialize screening 
    dshl = 0
    do iset=1,basis%nsets
      do ishell=1,basis%set(iset)%nshell
        dshl = dshl + 1
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
        shlrad(dshl) = sqrt(shell_gto_radius(alpha,d,sa%l,ntolnum))
      end do
    end do
  end subroutine

end module
