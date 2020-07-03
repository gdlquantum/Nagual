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
module nmintegrals
! Evaluation of gaussian AO integrals 
!
! Lit: S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!      M. Head-Gordon and J. A. Pople, J. Chem. Phys. 89, 5777 (1988)

  use nmparameter
  use nmtypes
  use nmfile
  use nmmatrix
  use nmmath
  use nmshell
  use nmbecke
  use nmunits

  implicit none

    public :: integrals_initialize
    public :: integrals_overlap
    public :: integrals_kinetic
    public :: integrals_core
    public :: integrals_ppi
    public :: integrals_ppi_three
    public :: integrals_ppi_two
    public :: integrals_ppi3
    public :: integrals_ppi2
    public :: integrals_dipole
    public :: integrals_set_ppi_type
    public :: non_zero_overlap
    public :: get_integration_pointers
    public :: psdaop

    private

      integer :: yuknqp
      parameter (yuknqp = 100)

      integer :: print_tape
      real(8) :: rho,coulw,yukw,yukg1,yukg2,yukv1,yukv2
      real(8) :: ftab(0:2*nmaxli+6,0:nmaxtgam)
      real(8) :: yukqp(yuknqp),yukqw(yuknqp)

contains

  subroutine integrals_ppi(ra,rb,rc,rd,sa,sb,sc,sd,ib,eps,skip,nsda,nsdb,nsdc)
  ! Evaluate electron Coulomb integrals
  ! Roberto Flores Moreno 2008, 2009, 2018
  implicit none
    integer :: nsda,nsdb,nsdc
    logical :: skip(2)
    real(8) :: eps
    real(8) :: ra(3),rb(3),rc(3),rd(3)
    real(8) :: ib(nmaxnco,nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc,sd

    integer :: ax,ay,az,bm1,bx,by,bz,cx,cy,cz,da,db,dc,dd,dleft,dright,dx,dy,dz
    integer :: i,ii,im1,im2,j,jj,jm1,k,l,lab,lcd,ll,lm1,lm2,lt,lx,ly,lz,m,mp1
    real(8) :: ab2,cd2,f1,f2,f3,f(0:2*nmaxli),zp,zp2,zq,zq2,zw,zw2
    real(8) :: xip,xiq,ep,eq,ew,pa,qc,wp,wq
    real(8) :: p(3),q(3),w(3),ecd(nmaxcon,nmaxcon)

    integer :: allocs
    real(8),allocatable :: v(:,:,:),hl(:,:),hr(:,:),h(:,:),hh(:,:,:)

    lab = sa%l + sb%l
    lcd = sc%l + sd%l
    lt = lab + lcd
    ab2 = sum((ra-rb)**2)
    cd2 = sum((rc-rd)**2)
    dleft = psdaop(lab)
    dright = psdaop(lcd)
    da = (sa%l+1)*(sa%l+2)/2
    db = (sb%l+1)*(sb%l+2)/2
    dc = (sc%l+1)*(sc%l+2)/2
    dd = (sd%l+1)*(sd%l+2)/2

    allocate(v(dleft,dright,0:lt),h(dleft,dright),stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: allocation failed')

    do k=1,sc%k
      do l=1,sd%k
        zq = sc%z(k) + sd%z(l)
        xiq = sc%z(k)*sd%z(l)/zq
        eq = xiq*cd2
        ecd(k,l) = exp(-eq)*sc%c(k)*sd%c(l)
      end do
    end do

    h(1:dleft,1:dright) = 0.0

    skip(1:2) = .true.
    do i=1,sa%k
      do j=1,sb%k
        zp = sa%z(i) + sb%z(j)
        zp2 = 2.0*zp
        xip = sa%z(i)*sb%z(j)/zp
        ep = xip*ab2
        p = (ra*(sa%z(i)) + rb*(sb%z(j)))/zp
        f1 = 2.0*exp(-ep)*sa%c(i)*sb%c(j)
        if (abs(2.0*pi*f1/zp).lt.eps) go to 2000
        f1 = f1*((sa%z(i))**nsda)*((sb%z(j))**nsdb)
        skip(1) = .false.
        do k=1,sc%k
          do l=1,sd%k
            zq = sc%z(k) + sd%z(l)
            zq2 = 2.0*zq
            xiq = sc%z(k)*sd%z(l)/zq
            eq = xiq*cd2
            q = (rc*(sc%z(k)) + rd*(sd%z(l)))/zq
            zw = zp + zq
            zw2 = 2*zw
            rho = zp*zq/zw
            ew = rho*sum((p-q)**2)
            w = (p*zp + q*zq)/zw
            if (abs(2.0*pi*ecd(k,l)/zq).lt.eps) go to 1000
            f2 = f1*ecd(k,l)*sqrt(rho/pi)*((pi**2/(zp*zq))**1.5)*((sc%z(k))**nsdc)
            skip(2) = .false.
            call gammaf(lt,ew,f)
            f(0:lt) = f2*f(0:lt)
            ! s-s-s-s type
            v(1,1,0:lt) = f(0:lt)
            ! vrr right side
            do ll=1,lcd
              lm1 = ll - 1
              lm2 = max(0,ll-2)
              m = lt - ll
              mp1 = m + 1
              qc = q(3) - rc(3)
              wq = w(3) - q(3)
              ii = gaop(0,0,ll)
              im1 = gaop(0,0,lm1)
              im2 = gaop(0,0,lm2)
              v(1,ii,0:m) = qc*v(1,im1,0:m) + wq*v(1,im1,1:mp1) +       &
     &        (float(lm1)/zq2)*(v(1,im2,0:m)-zp/zw*v(1,im2,1:mp1))
              qc = q(2) - rc(2)
              wq = w(2) - q(2)
              do ly=1,ll
                lm1 = ly - 1
                lm2 = max(0,ly-2)
                lz = ll - ly
                ii = gaop(0,ly,lz)
                im1 = gaop(0,lm1,lz)
                im2 = gaop(0,lm2,lz)
                v(1,ii,0:m) = qc*v(1,im1,0:m) + wq*v(1,im1,1:mp1) +     &
     &          (float(lm1)/zq2)*(v(1,im2,0:m)-zp/zw*v(1,im2,1:mp1))
              end do
              qc = q(1) - rc(1)
              wq = w(1) - q(1)
              do lx=1,ll
                lm1 = lx - 1
                lm2 = max(0,lx-2)
                do ly=0,ll-lx
                  lz = ll - lx - ly
                  ii = gaop(lx,ly,lz)
                  im1 = gaop(lm1,ly,lz)
                  im2 = gaop(lm2,ly,lz)
                  v(1,ii,0:m) = qc*v(1,im1,0:m) + wq*v(1,im1,1:mp1) +   &
     &            (float(lm1)/zq2)*(v(1,im2,0:m)-zp/zw*v(1,im2,1:mp1))
                end do
              end do
            end do
           ! vrr left side
            do ll=1,lab
              lm1 = ll - 1
              lm2 = max(0,ll-2)
              m = lab - ll
              mp1 = m + 1
              pa = p(3)-ra(3)
              wp = w(3)-p(3)
              ii = gaop(0,0,ll)
              im1 = gaop(0,0,lm1)
              im2 = gaop(0,0,lm2)
              v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*         &
     &        v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*(                &
     &        v(im2,1:dright,0:m)-zq/zw*v(im2,1:dright,1:mp1))
              do bz=1,lcd
                f3 = float(bz)/zw2
                bm1 = bz - 1
                do bx=0,lcd-bz
                  do by=0,lcd-bx-bz
                    jj = gaop(bx,by,bz)
                    jm1 = gaop(bx,by,bm1)
                    v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                  end do
                end do
              end do
              pa = p(2)-ra(2)
              wp = w(2)-p(2)
              do ly=1,ll
                lm1 = ly - 1
                lm2 = max(0,ly-2)
                lz = ll - ly
                ii = gaop(0,ly,lz)
                im1 = gaop(0,lm1,lz)
                im2 = gaop(0,lm2,lz)
                v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*       &
     &          v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*(              &
     &          v(im2,1:dright,0:m)-zq/zw*v(im2,1:dright,1:mp1))
                do by=1,lcd
                  f3 = float(by)/zw2
                  bm1 = by - 1
                  do bx=0,lcd-by
                    do bz=0,lcd-bx-by
                      jj = gaop(bx,by,bz)
                      jm1 = gaop(bx,bm1,bz)
                      v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                    end do
                  end do
                end do
              end do
              pa = p(1)-ra(1)
              wp = w(1)-p(1)
              do lx=1,ll
                lm1 = lx - 1
                lm2 = max(0,lx-2)
                do ly=0,ll-lx
                  lz = ll - lx - ly
                  ii = gaop(lx,ly,lz)
                  im1 = gaop(lm1,ly,lz)
                  im2 = gaop(lm2,ly,lz)
                  v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*     &
     &            v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*(            &
     &            v(im2,1:dright,0:m)-zq/zw*v(im2,1:dright,1:mp1))
                  do bx=1,lcd
                    f3 = float(bx)/zw2
                    bm1 = bx - 1
                    do by=0,lcd-bx
                      do bz=0,lcd-bx-by
                        jj = gaop(bx,by,bz)
                        jm1 = gaop(bm1,by,bz)
                        v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                      end do
                    end do
                  end do
                end do
              end do
            end do
           ! Contraction
            h(1:dleft,1:dright) = h(1:dleft,1:dright) + v(1:dleft,1:dright,0)
 1000       continue
          end do
        end do
 2000       continue
      end do
    end do
    deallocate(v,stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: deallocation failed')
    if (skip(2)) then
      deallocate(h,stat=allocs)
      if (allocs.gt.0) call file_error('integrals_ppi: deallocation failed')
      ib(:da,:db,:dc,:dd) = 0.0
      return
    end if
    allocate(hr(dright,dright),hh(dleft,dc,dd),stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: allocation failed')
    do i=1,dleft
      hr(1:dright,1) = h(i,1:dright)
      call hrr(sc%l,sd%l,rd-rc,hr,dright)
      do cx=0,sc%l
        do cy=0,sc%l-cx
          cz = sc%l - cx - cy
          do dx=0,sd%l
            do dy=0,sd%l-dx
              dz = sd%l-dx-dy
              hh(i,raop(cx,cy,cz),raop(dx,dy,dz)) =                     &
     &        hr(gaop(cx,cy,cz),gaop(dx,dy,dz))
            end do
         end do
        end do
      end do
    end do
    deallocate(h,hr,stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: deallocation failed')
    allocate(hl(dleft,dleft),stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: allocation failed')
    do cx=0,sc%l
      do cy=0,sc%l-cx
        cz = sc%l - cx - cy
        k = raop(cx,cy,cz)
        do dx=0,sd%l
          do dy=0,sd%l-dx
            dz = sd%l-dx-dy
            l = raop(dx,dy,dz)
            hl(1:dleft,1) = hh(1:dleft,k,l)
            call hrr(sa%l,sb%l,rb-ra,hl,dleft)
            do ax=0,sa%l
              do ay=0,sa%l-ax
                az = sa%l - ax - ay
                i = raop(ax,ay,az)
                do bx=0,sb%l
                  do by=0,sb%l-bx
                    bz = sb%l-bx-by
                    j = raop(bx,by,bz)
                    ib(i,j,k,l) = hl(gaop(ax,ay,az),gaop(bx,by,bz))*    &
     &              sa%norm(i)*sb%norm(j)*sc%norm(k)*sd%norm(l)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(hl,hh,stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi: deallocation failed')
  end subroutine 

  subroutine integrals_core(ra,rb,rc,sa,sb,ib,nsda,nsdb)
  ! Evaluate nuclear atraction integrals
  ! Roberto Flores Moreno 2008
  implicit none
    integer :: nsda,nsdb
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: ib(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb

    integer :: ax,ay,az,bx,by,bz,daop,i,j,k,l,lab,mm,ma
    real(8) :: factor,zp,xi,sf,r2,t
    real(8) :: r(3),p(3),c(3),q(3),f(0:nmaxli)

    integer :: allocs
    real(8),allocatable :: v(:,:),h(:,:)

    lab = sa%l + sb%l
    daop = psdaop(lab)

    allocate(v(daop,0:lab),h(daop,daop),stat=allocs)
    if (allocs.gt.0) call file_error('core_integrals: allocation failed')

    ! Initialize
    r = rb - ra
    r2 = sum(r**2)

    h(1:daop,1) = 0.0
    do i=1,sa%k
      do j=1,sb%k
        zp = sa%z(i) + sb%z(j)
        xi = sa%z(i)*sb%z(j)/zp
        sf = sb%z(j)/zp
        p = (ra*(sa%z(i)) + rb*(sb%z(j)))/zp
        factor = 2.0*pi/zp*exp(-xi*r2)*sa%c(i)*sb%c(j)
        q = p - rc
        t = zp*sum(q**2)
        ! s-s type
        call gamma_coulomb(lab,t,f)
        do k=0,lab
          v(gaop(0,0,0),k) = factor*f(k)
        end do
        ! ax up
        mm = lab
        if (mm.gt.0) then
          v(gaop(1,0,0),0:mm-1) = sf*r(1)*v(gaop(0,0,0),0:mm-1)         &
     &                          - q(1)*v(gaop(0,0,0),1:mm)
        end if
        do ax=2,mm
          ma = mm - ax
          v(gaop(ax,0,0),0:ma) = sf*r(1)*v(gaop(ax-1,0,0),0:ma)         &
     &    -q(1)*v(gaop(ax-1,0,0),1:ma+1) + (ax-1)/(2.0*zp)*             &
     &    (v(gaop(ax-2,0,0),0:ma)-v(gaop(ax-2,0,0),1:ma+1))
        end do
        ! ay up
        do ax=0,lab
          mm = lab - ax
          if (mm.gt.0) then
            v(gaop(ax,1,0),0:mm-1) = sf*r(2)*v(gaop(ax,0,0),0:mm-1)     &
     &                             - q(2)*v(gaop(ax,0,0),1:mm)
          end if
          do ay=2,mm
            ma = mm - ay
            v(gaop(ax,ay,0),0:ma) = sf*r(2)*v(gaop(ax,ay-1,0),0:ma)     &
     &      -q(2)*v(gaop(ax,ay-1,0),1:ma+1) + (ay-1)/(2.0*zp)*          &
     &      (v(gaop(ax,ay-2,0),0:ma)-v(gaop(ax,ay-2,0),1:ma+1))
          end do
        end do
        ! az up
        do ax=0,lab
          do ay=0,lab-ax
            mm = lab - ax - ay
            if (mm.gt.0) then
              v(gaop(ax,ay,1),0:mm-1) = sf*r(3)*v(gaop(ax,ay,0),0:mm-1) &
     &                                - q(3)*v(gaop(ax,ay,0),1:mm)
            end if
            do az=2,mm
              ma = mm - az
              v(gaop(ax,ay,az),0:ma) = sf*r(3)*v(gaop(ax,ay,az-1),0:ma) &
     &        -q(3)*v(gaop(ax,ay,az-1),1:ma+1) + (az-1)/(2.0*zp)*       &
     &        (v(gaop(ax,ay,az-2),0:ma)-v(gaop(ax,ay,az-2),1:ma+1))
            end do
          end do
        end do
        ! contraction
        h(1:daop,1) = h(1:daop,1) + v(1:daop,0)*(sa%z(i)**nsda)*(sb%z(j)**nsdb)
      end do
    end do
    ! horizontal recurrence relation
    call hrr(sa%l,sb%l,r,h,daop)
    ! save and normalize
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        do bx=0,sb%l
          do by=0,sb%l-bx
            bz = sb%l-bx-by
            ib(raop(ax,ay,az),raop(bx,by,bz)) = sa%norm(raop(ax,ay,az))*&
     &      sb%norm(raop(bx,by,bz))*h(gaop(ax,ay,az),gaop(bx,by,bz))
          end do
        end do
      end do
    end do
    deallocate(v,h,stat=allocs)
    if (allocs.gt.0) call file_error('core_integrals: deallocation failed')
  end subroutine 

  subroutine integrals_kinetic(ra,rb,sa,sb,ib,nsd)
  ! Evaluate kinetic energy integrals
  ! Roberto Flores Moreno,  2008, 2014, 2018
  implicit none
    integer :: nsd
    real(8) :: ra(3),rb(3),ib(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb

    integer :: ap,ax,ay,az,bx,by,bz,na,nb
    type(nshell) :: wsa,wsb

    integer :: allocs,dover
    real(8),allocatable :: ob(:,:)

    dover = ((sa%l+sb%l+nsd+3)*(sa%l+sb%l+nsd+4))/2
    allocate(ob(dover,dover),stat=allocs)
    if (allocs.gt.0) call file_error('kinetic_integrals: allocation failed')

    na = (sa%l+1)*(sa%l+2)/2
    nb = (sb%l+1)*(sb%l+2)/2
    wsa = sa
    wsb = sb
    wsa%norm(:) = 1.0
    wsb%norm(:) = 1.0

    call integrals_overlap(ra,rb,wsa,wsb,ob,dover,1+nsd)
    ib(1:na,1:nb) = (2*sa%l+3)*ob(1:na,1:nb)
    wsa%l = sa%l + 2
    call integrals_overlap(ra,rb,wsa,wsb,ob,dover,2+nsd)
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        ap = raop(ax,ay,az)
        ib(ap,1:nb) = ib(ap,1:nb) - 2.0*ob(raop(ax+2,ay,az),1:nb)       &
     &                            - 2.0*ob(raop(ax,ay+2,az),1:nb)       &
     &                            - 2.0*ob(raop(ax,ay,az+2),1:nb)              
      end do
    end do
    if (sa%l.gt.1) then
      wsa%l = sa%l - 2
      call integrals_overlap(ra,rb,wsa,wsb,ob,dover,0+nsd)
      do ax=0,sa%l
        do ay=0,sa%l-ax
          az = sa%l - ax - ay
          ap = raop(ax,ay,az)
          if (ax.gt.1) then
            ib(ap,1:nb) = ib(ap,1:nb)-(ax*(ax-1)/2.0)*ob(raop(ax-2,ay,az),1:nb)
          end if
          if (ay.gt.1) then
            ib(ap,1:nb) = ib(ap,1:nb)-(ay*(ay-1)/2.0)*ob(raop(ax,ay-2,az),1:nb) 
          end if
          if (az.gt.1) then
            ib(ap,1:nb) = ib(ap,1:nb)-(az*(az-1)/2.0)*ob(raop(ax,ay,az-2),1:nb) 
          end if
        end do
      end do
    end if

    deallocate(ob,stat=allocs)
    if (allocs.gt.0) call file_error('kinetic_integrals: deallocation failed')

    ! normalize
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        do bx=0,sb%l
          do by=0,sb%l-bx
            bz = sb%l-bx-by
            ib(raop(ax,ay,az),raop(bx,by,bz)) = sa%norm(raop(ax,ay,az))*&
     &      sb%norm(raop(bx,by,bz))*ib(raop(ax,ay,az),raop(bx,by,bz))
          end do
        end do
      end do
    end do
  end subroutine 

  subroutine integrals_overlap(ra,rb,sa,sb,ib,db,nsd)
  ! Evaluate overlap integrals
  ! Roberto Flores-Moreno, 2008
  implicit none
    integer :: db,nsd
    real(8) :: ra(3),rb(3),ib(db,db)
    type(nshell) :: sa,sb

    integer :: daop,i,j,l,lab,ax,ay,az,bx,by,bz
    real(8) :: zp,zp2,xi,r2
    real(8) :: r(3),rr(3)

    integer :: allocs
    real(8),allocatable :: v(:),h(:,:)

    lab = sa%l + sb%l
    daop = psdaop(lab)

    allocate(v(daop),h(daop,daop),stat=allocs)
    if (allocs.gt.0) call file_error('overlap_integrals: allocation failed')

    ! Initialize
    r = rb-ra
    r2 = sum(r**2)

    h(1:daop,1) = 0.0

    ! Main loop on primitives
    do i=1,sa%k
      do j=1,sb%k
        zp = sa%z(i) + sb%z(j)
        zp2 = 2.0*zp
        xi = sa%z(i)*sb%z(j)/zp
        rr = sb%z(j)/zp*r
        ! s-s type
        v(gaop(0,0,0)) = ((pi/zp)**1.5)*exp(-xi*r2)
        ! vertical recurrence relation x
        if (lab.gt.0) v(gaop(1,0,0)) = rr(1)*v(gaop(0,0,0))
        do ax=2,lab
          v(gaop(ax,0,0)) = rr(1)*v(gaop(ax-1,0,0)) +                   &
     &                (float(ax-1)/zp2)*v(gaop(ax-2,0,0))
        end do
        ! vertical recurrence relation y
        do ax=0,lab
          if (lab-ax.gt.0) v(gaop(ax,1,0)) = rr(2)*v(gaop(ax,0,0))
          do ay=2,lab-ax
            v(gaop(ax,ay,0)) = rr(2)*v(gaop(ax,ay-1,0)) +               &
     &                        ((ay-1)/zp2)*v(gaop(ax,ay-2,0))
          end do
        end do
        ! vertical recurrence relation z
        do ax=0,lab
          do ay=0,lab-ax
            if (lab-ax-ay.gt.0) v(gaop(ax,ay,1)) = rr(3)*v(gaop(ax,ay,0))
            do az=2,lab-ax-ay
              v(gaop(ax,ay,az)) = rr(3)*v(gaop(ax,ay,az-1)) +           &
     &                           ((az-1)/zp2)*v(gaop(ax,ay,az-2))
            end do
          end do
        end do
        ! contraction
        h(1:daop,1) = h(1:daop,1) + (sa%z(i)**nsd)*sa%c(i)*sb%c(j)*v(1:daop)
      end do
    end do
    ! horizontal recurrence relation
    call hrr(sa%l,sb%l,r,h,daop)
    ! save and normalize
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        do bx=0,sb%l
          do by=0,sb%l-bx
            bz = sb%l-bx-by
            ib(raop(ax,ay,az),raop(bx,by,bz)) = sa%norm(raop(ax,ay,az))*&
     &      sb%norm(raop(bx,by,bz))*h(gaop(ax,ay,az),gaop(bx,by,bz))
          end do
        end do
      end do
    end do
    
    deallocate(v,h,stat=allocs)
    if (allocs.gt.0) call file_error('overlap_integrals: deallocation failed')
  end subroutine 

  subroutine hrr(la,lb,r,h,dh)
  ! Calculation of integrals over horizontal recurrence relation.
  ! Roberto Flores-Moreno, 2008
  implicit none
    integer :: dh,la,lb
    real(8) :: r(3),h(dh,dh)

    integer :: ax,ay,az,bx,by,bz,lab,l1,l2

    lab = la + lb
    do l1=1,lb
      do l2=la,lab-l1
        do ax=0,l2
          do ay=0,l2-ax
            az = l2 - ax - ay
            ! bz up
            h(gaop(ax,ay,az),gaop(0,0,l1)) =                            &
     &      h(gaop(ax,ay,az+1),gaop(0,0,l1-1)) -                        &
     &      r(3)*h(gaop(ax,ay,az),gaop(0,0,l1-1))
            ! by up
            do by=1,l1
              bz = l1 - by
              h(gaop(ax,ay,az),gaop(0,by,bz)) =                         &
     &        h(gaop(ax,ay+1,az),gaop(0,by-1,bz)) -                     &
     &        r(2)*h(gaop(ax,ay,az),gaop(0,by-1,bz))
            end do
            ! bx up
            do bx=1,l1
              do by=0,l1-bx
                bz = l1 - bx - by
                h(gaop(ax,ay,az),gaop(bx,by,bz)) =                      &
     &          h(gaop(ax+1,ay,az),gaop(bx-1,by,bz)) -                  &
     &          r(1)*h(gaop(ax,ay,az),gaop(bx-1,by,bz))
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine hrr

  subroutine integrals_initialize(m,tape)
  ! Initialize integration routines
  ! Roberto Flores-Moreno, 2008
  implicit none
    type(nmolecule) :: m
    integer :: tape

    integer :: i,iatom,ishell,j,k,l,lx,ly,lz
    real(8) :: factor,norm
    type(nshell) :: s

    integer :: nitermax,nmax
    real(8) :: eps
    real(8) :: bessel,expterm,prefak,preterm,produkt,serie,sumterm,term,ttab

    integer :: allocs
    real(8),allocatable :: r(:)

    print_tape = tape

    call shell_initialize_pointers

    ! Gamma function
    nitermax = 30 
    nmax = 2*nmaxli+6
    eps = 1.0e-15
    allocate(r(nitermax+10),stat=allocs)
    do i=0,nmax
      ftab(i,0) = 1.0/(2*i+1)
    end do
    do i=1,nmaxtgam
      ttab = float(i)/10.0
      r(nitermax+10) = 0.0
      do j=1,nitermax+9
        r(nitermax+10-j) = - ttab/(4*(nitermax+10-j) + 2.0 - ttab*r(nitermax+11-j))
      end do
      bessel = (2*sinh(ttab/2))/ttab
      prefak = exp(-ttab/2)*bessel
      term = 1.0
      serie = prefak*(1.0/(2.0*nmax + 1.0))
      do k=1,nitermax
        preterm = (2.0*k + 1.0)/(2.0*nmax + 1.0)
        term = term*(2.0*nmax - 2.0*k + 1.0)/(2.0*nmax + 2.0*k + 1.0)
        produkt = 1.0
        do l=1,k
          produkt = produkt*r(l)
        end do
        sumterm = prefak*preterm*term*produkt
        if (abs(sumterm).le.eps) then
          go to 500
        else
          serie = serie + sumterm
        end if
      end do
      call file_error('fatal integrals_initialize, contact support')
  500 continue
      ftab(nmax,i) = serie
      expterm = exp(-ttab)
      do j=1,nmax
        ftab(nmax-j,i) = 1.0/(2*(nmax-j)+1)*(2*ttab*ftab(nmax+1-j,i) + expterm)
      end do
    end do
    deallocate(r,stat=allocs)
  end subroutine 

  subroutine integrals_set_ppi_type(simple)
  ! Set particle-particle interaction type
  ! Roberto Flores-Moreno, 2018
  implicit none
    real(8) :: simple(*)

    coulw = simple(1)
    yukw = simple(2)

    if (yukw.ne.0.0) then
      yukv1 = units_ev_to_au(2.0e+9)
      yukg1 = 1.0/units_angstrom_to_bohr(0.2e-5)
      yukv2 = units_ev_to_au(9.525e+7)
      yukg2 = 1.0/units_angstrom_to_bohr(1.68e-5)
      call becke_tskgc(yukqp,yukqw,yuknqp)
    end if
  end subroutine

  subroutine gammaf(m,t,f)
  ! Calculation of Boys and related functions
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: m
    real(8) :: t
    real(8) :: f(0:*)

    real(8) :: g(0:2*nmaxli)

    f(0:m) = 0.0

    ! Coulombic interaction
    if (coulw.ne.0.0) then
      call gamma_coulomb(m,t,g)
      f(0:m) = f(0:m) + coulw*g(0:m)
    end if

    if (yukw.ne.0.0) then
      call gamma_yukawa(m,t,g)
      f(0:m) = f(0:m) + yukw*g(0:m)
    end if
  end subroutine

  subroutine gamma_coulomb(m,t,f)
  ! Calculation of the incomplete gamma function F(t)
  ! for ppi integrals over Gaussian functions.
  ! L. E. McMurchie and E.R. Davidson, J. Comp. Phys. 26, 218 (1978).
  ! Original version from Andreas M. Koester, 1996
  ! Roberto Flores-Moreno, 2018
  implicit none
    real(8) :: eps

    integer :: m
    real(8) :: t
    real(8) :: f(0:*)

    integer :: i,k,ttab
    real(8) :: a,b,c,d,expterm

    eps = 1.0e-13
    if (t.lt.0.0) t = eps
    if (t.le.eps) then
      f(m) = 1.0/(2.0*m + 1.0)
      do i=1,m
        f(m-i) = 1.0/(2.0*(m-i) + 1.0)
      end do
      return
    else if (t.le.12.0) then
      ttab = nint(10*t)
      f(m) = ftab(m,ttab)
      do k=1,6
        f(m) = f(m) + ftab(m+k,ttab)*(float(ttab)/10.0 - t)**k/math_factorial(k)
      end do
      if (m.gt.0) expterm = exp(-t)
      do i=1,m
        f(m-i) = 1.0/(2*(m-i)+1) * (2*t*f(m+1-i) + expterm)
      end do
      return
    else if (t.le.15.0) then
      a = 0.4999489092
      b = 0.2473631686
      c = 0.3211809090
      d = 0.3811559346
      f(0) = 0.5*sqrt(pi/t) - (exp(-t)/t)*(a - b/t + c/t**2 - d/t**3)
    else if (t.le.18.0) then
      a = 0.4998436875
      b = 0.2424943800
      c = 0.2464284500
      f(0) = 0.5*sqrt(pi/t) - (exp(-t)/t)*(a - b/t + c/t**2)
    else if (t.le.24.0) then
      a = 0.4990931620
      b = 0.2152832000
      f(0) = 0.5*sqrt(pi/t) - (exp(-t)/t)*(a - b/t)
    else if (t.le.30.0) then
      a = 0.49000000
      f(0) = 0.5*sqrt(pi/t) - (exp(-t)/t)*a
    else
      f(0) = 0.5*sqrt(pi/t)
    end if
    if (t.gt.(2.0*m+36)) then
      do i=1,m
        f(i) = (2*i-1)/(2*t)*f(i-1)
      end do
    else
      expterm = exp(-t)
      do i=1,m
        f(i) = 1/(2*t) * ((2*i-1)*f(i-1) - expterm)
      end do
    end if
  end subroutine 

  subroutine gamma_yukawa(m,t,g)
  ! Calculation of extended Boys function
  ! Henry Gonzalez & R. Flores-Moreno, 2018
  implicit none
    integer :: m
    real(8) :: t
    real(8) :: g(0:*)

    integer :: i,iqp
    real(8) :: a1,a2,factor,factor0,expo1,expo2,y,y2,wy,u1,u2

    u1 = 0.25*yukg1*yukg1/rho
    u2 = 0.25*yukg2*yukg2/rho
    a1 = (yukv1/yukg1)!*exp(u1)
    a2 = (yukv2/yukg2)!*exp(u2)

    g(0:m) = 0.0
    do iqp=1,yuknqp
      y = 0.5*yukqp(iqp) + 0.5!float(iqp)/float(nqp+1)
      wy = 0.5*yukqw(iqp)!1.0/float(nqp)
      y2 = y*y
      expo1 = - u1*(1.0 - 1.0/y2) + t*y2
      expo2 = - u2*(1.0 - 1.0/y2) + t*y2
      factor = 0.0
      if (abs(expo1).lt.1.0e-7) then
        factor = factor + a1*(1.0-expo1)
      else
        factor = factor + a1*exp(-expo1)
      end if
      if (abs(expo2).lt.1.0e-7) then
        factor = factor - a2*(1.0-expo2)
      else
        factor = factor - a2*exp(-expo2)
      end if
      factor = wy*factor
      g(0) = g(0) + factor
      do i=1,m
        factor = factor*y2
        g(i) = g(i) + factor
      end do
    end do
  end subroutine

  integer function psdaop(l)
  ! Get physical storage dimension for aop
  ! Roberto Flores-Moreno, 2008
  implicit none
    integer :: l

    psdaop = (l**3-l)/6 + l**2 + 2*l + 1
  end function psdaop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Screening

  logical function non_zero_overlap(ra,rb,sa,sb,eps)
  ! Determine if there is non-zero overlap among two shells
  ! Roberto Flores-Moreno, Aug 2008
  implicit none
    real(8) :: eps
    real(8) :: ra(3),rb(3)
    type(nshell) :: sa,sb

    integer :: i,j
    real(8) :: ab2,ep,gp,xi,zp

    non_zero_overlap = .false.
    ab2 = sum((ra-rb)**2)

    do i=1,sa%k
      do j=1,sb%k
        zp = sa%z(i) + sb%z(j)
        xi = sa%z(i)*sb%z(j)/zp
        ep = xi*ab2
        gp = 2.0*exp(-ep)*sa%c(i)*sb%c(j)
        if (abs(2.0*pi*gp/zp).ge.eps) then
          non_zero_overlap = .true.
          return
        end if
      end do
    end do
  end function

  subroutine get_integration_pointers(rraop,rgaop)
  ! Get integration pointers for compatibility of other modules
  ! Roberto Flores-Moreno, Oct 2008
  implicit none
    integer :: rgaop(0:nmaxli,0:nmaxli,0:nmaxli)
    integer :: rraop(0:nmaxli,0:nmaxli,0:nmaxli)
    
    rgaop(0:nmaxli,0:nmaxli,0:nmaxli) = gaop(0:nmaxli,0:nmaxli,0:nmaxli)
    rraop(0:nmaxli,0:nmaxli,0:nmaxli) = raop(0:nmaxli,0:nmaxli,0:nmaxli)
  end subroutine

  subroutine integrals_ppi_three(ra,rb,rc,sa,sb,sc,ib,tol,skip,nsda,nsdb)
  ! Evaluate three-centerelectron Coulomb integrals
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    integer :: nsda,nsdb
    logical :: skip
    real(8) :: tol
    real(8) :: ra(3),rb(3),rc(3)
    type(nshell) :: sa,sb,sc
    real(8) :: ib(nmaxnco,nmaxnco,nmaxnco)

    integer :: na,nb,nc
    logical :: skipd(2)
    real(8) :: rd(3)
    type(nshell) :: sd
    real(8) :: ob(nmaxnco,nmaxnco,nmaxnco,nmaxnco)

    rd = rc
    sd%l = 0
    sd%norm = 1.0
    sd%k = 1
    sd%c(1) = 1.0
    sd%z(1) = 0.0

    call integrals_ppi(ra,rb,rc,rd,sa,sb,sc,sd,ob,tol,skipd,nsda,nsdb,0)

    skip = skipd(1)
    na = (sa%l+1)*(sa%l+2)/2
    nb = (sb%l+1)*(sb%l+2)/2
    nc = (sc%l+1)*(sc%l+2)/2
    ib(1:na,1:nb,1:nc) = ob(1:na,1:nb,1:nc,1)
  end subroutine

  subroutine integrals_ppi_two(ra,rc,sa,sc,ib,tol,nsd)
  ! Evaluate tw-center electron Coulomb integrals
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    integer :: nsd
    real(8) :: tol
    real(8) :: ra(3),rc(3)
    type(nshell) :: sa,sc
    real(8) :: ib(nmaxnco,nmaxnco)

    integer :: na,nc
    logical :: skip(2)
    real(8) :: rb(3),rd(3)
    type(nshell) :: sb,sd
    real(8) :: ob(nmaxnco,nmaxnco,nmaxnco,nmaxnco)

    rb = ra
    sb%l = 0
    sb%norm = 1.0
    sb%k = 1
    sb%c(1) = 1.0
    sb%z(1) = 0.0

    rd = rc
    sd%l = 0
    sd%norm = 1.0
    sd%k = 1
    sd%c(1) = 1.0
    sd%z(1) = 0.0

    call integrals_ppi(ra,rb,rc,rd,sa,sb,sc,sd,ob,tol,skip,nsd,0,0)

    na = (sa%l+1)*(sa%l+2)/2
    nc = (sc%l+1)*(sc%l+2)/2
    ib(1:na,1:nc) = ob(1:na,1,1:nc,1)
  end subroutine

  subroutine integrals_ppi3(ra,rb,rc,sa,sb,sc,ib,eps,skip,nsda,nsdb,&
                            h,hl,v,left,right,ltmax)
  ! Evaluate 3-index particle-particle interaction integrals
  ! Roberto Flores Moreno, 2018
  implicit none
    integer :: left,ltmax,nsda,nsdb,right
    logical :: skip
    real(8) :: eps
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: ib(nmaxnco,nmaxnco,nmaxnco)
    real(8) :: h(left,right),hl(left,left),v(left,right,0:ltmax)
    type(nshell) :: sa,sb,sc

    integer :: ax,ay,az,ba,bb,bc,bm1,bx,by,bz,cx,cy,cz,da,db,dc,dleft,dright
    integer :: i,ii,im1,im2,j,jj,jm1,k,kk,lab,lc,ll,lm1,lm2,lt,lx,ly,lz,m,mp1
    real(8) :: ab2,cd2,fz,f0,f1,f2,f3,f(0:2*nmaxli),zp,zp2,zq,zq2,zw,zw2
    real(8) :: xip,xiq,ep,eq,ew,pa,qc,wp,wq
    real(8) :: p(3),q(3),w(3)

    lab = sa%l + sb%l
    lc = sc%l 
    lt = lab + lc
    ab2 = sum((ra-rb)**2)
    dleft = psdaop(lab)
    dright = psdaop(lc)
    da = (sa%l+1)*(sa%l+2)/2
    db = (sb%l+1)*(sb%l+2)/2
    dc = (sc%l+1)*(sc%l+2)/2

    h(1:dleft,1:dright) = 0.0

    skip = .true.
    do i=1,sa%k
      do j=1,sb%k
        zp = sa%z(i) + sb%z(j)
        zp2 = zp + zp
        xip = sa%z(i)*sb%z(j)/zp
        ep = xip*ab2
        p = (ra*(sa%z(i)) + rb*(sb%z(j)))/zp
        f1 = 2.0*exp(-ep)*sa%c(i)*sb%c(j)
        if (abs(2.0*pi*f1/zp).lt.eps) go to 2000
        f1 = f1*((sa%z(i))**nsda)*((sb%z(j))**nsdb)
        skip = .false.
        !rfm do k=1,sc%k
          zq = sc%z(1) 
          zq2 = zq + zq
          q = rc
          zw = zp + zq
          zw2 = zw + zw
          rho = zp*zq/zw
          ew = rho*sum((p-q)**2)
          w = (p*zp + q*zq)/zw
          f2 = f1*sc%c(1)*sqrt(rho/pi)*((pi**2/(zp*zq))**1.5)
          call gammaf(lt,ew,f)
          f(0:lt) = f2*f(0:lt)
          ! s-s-s-s type
          v(1,1,0:lt) = f(0:lt)
          ! vrr right side
          fz = zp/zw
          do ll=1,lc
            lm1 = ll - 1
            m = lt - ll
            mp1 = m + 1
            wq = w(3) - q(3)
            ii = gaop(0,0,ll)
            im1 = gaop(0,0,lm1)
            if (lm1.eq.0) then
              v(1,ii,0:m) = wq*v(1,im1,1:mp1)
            else
              lm2 = ll-2
              im2 = gaop(0,0,lm2)
              v(1,ii,0:m) = wq*v(1,im1,1:mp1) +       &
     &        (float(lm1)/zq2)*(v(1,im2,0:m)-fz*v(1,im2,1:mp1))
            end if
            wq = w(2) - q(2)
            do ly=1,ll
              lm1 = ly - 1
              lz = ll - ly
              ii = gaop(0,ly,lz)
              im1 = gaop(0,lm1,lz)
              if (lm1.eq.0) then
                v(1,ii,0:m) = wq*v(1,im1,1:mp1)
              else
                lm2 = ly-2
                im2 = gaop(0,lm2,lz)
                v(1,ii,0:m) = wq*v(1,im1,1:mp1) +     &
     &          (float(lm1)/zq2)*(v(1,im2,0:m)-fz*v(1,im2,1:mp1))
              end if
            end do
            wq = w(1) - q(1)
            do lx=1,ll
              lm1 = lx - 1
              if (lm1.gt.0) f0 = float(lm1)/zq2
              lm2 = max(0,lx-2)
              do ly=0,ll-lx
                lz = ll - lx - ly
                ii = gaop(lx,ly,lz)
                im1 = gaop(lm1,ly,lz)
                if (lm1.eq.0) then
                  v(1,ii,0:m) = wq*v(1,im1,1:mp1) 
                else
                  im2 = gaop(lm2,ly,lz)
                  v(1,ii,0:m) = wq*v(1,im1,1:mp1) +   &
     &            f0*(v(1,im2,0:m)-fz*v(1,im2,1:mp1))
                end if
              end do
            end do
          end do
         ! vrr left side
          fz = zq/zw
          do ll=1,lab
            lm1 = ll - 1
            m = lab - ll
            mp1 = m + 1
            pa = p(3)-ra(3)
            wp = w(3)-p(3)
            ii = gaop(0,0,ll)
            im1 = gaop(0,0,lm1)
            if (lm1.eq.0) then
              v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*         &
     &        v(im1,1:dright,1:mp1) 
            else
              lm2 = ll-2
              im2 = gaop(0,0,lm2)
              v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*         &
     &        v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*(                &
     &        v(im2,1:dright,0:m)-fz*v(im2,1:dright,1:mp1))
            end if
            do bz=1,lc
              f3 = float(bz)/zw2
              bm1 = bz - 1
              do bx=0,lc-bz
                do by=0,lc-bx-bz
                  jj = gaop(bx,by,bz)
                  jm1 = gaop(bx,by,bm1)
                  v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                end do
              end do
            end do
            pa = p(2)-ra(2)
            wp = w(2)-p(2)
            do ly=1,ll
              lm1 = ly - 1
              if (lm1.gt.0) f0 = float(lm1)/zp2
              lz = ll - ly
              ii = gaop(0,ly,lz)
              im1 = gaop(0,lm1,lz)
              if (lm1.eq.0) then
                v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*       &
     &          v(im1,1:dright,1:mp1) 
              else
                lm2 = ly-2
                im2 = gaop(0,lm2,lz)
                v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*       &
     &          v(im1,1:dright,1:mp1) + f0*(              &
     &          v(im2,1:dright,0:m)-fz*v(im2,1:dright,1:mp1))
              end if
              do by=1,lc
                f3 = float(by)/zw2
                bm1 = by - 1
                do bx=0,lc-by
                  do bz=0,lc-bx-by
                    jj = gaop(bx,by,bz)
                    jm1 = gaop(bx,bm1,bz)
                    v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                  end do
                end do
              end do
            end do
            pa = p(1)-ra(1)
            wp = w(1)-p(1)
            do lx=1,ll
              lm1 = lx - 1
              if (lm1.gt.0) f0 = float(lm1)/zp2
              lm2 = max(0,lx-2)
              do ly=0,ll-lx
                lz = ll - lx - ly
                ii = gaop(lx,ly,lz)
                im1 = gaop(lm1,ly,lz)
                if (lm1.eq.0) then
                  v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*     &
     &            v(im1,1:dright,1:mp1) 
                else
                  im2 = gaop(lm2,ly,lz)
                  v(ii,1:dright,0:m) = pa*v(im1,1:dright,0:m) + wp*     &
     &            v(im1,1:dright,1:mp1) + f0*(            &
     &            v(im2,1:dright,0:m)-fz*v(im2,1:dright,1:mp1))
                end if
                do bx=1,lc
                  f3 = float(bx)/zw2
                  bm1 = bx - 1
                  do by=0,lc-bx
                    do bz=0,lc-bx-by
                      jj = gaop(bx,by,bz)
                      jm1 = gaop(bm1,by,bz)
                      v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
                    end do
                  end do
                end do
              end do
            end do
          end do
          ! Contraction
          h(1:dleft,1:dright) = h(1:dleft,1:dright) + v(1:dleft,1:dright,0)
 1000     continue
        !rfm end do
 2000   continue
      end do
    end do
    q = rb - ra
    ba = shell_dimgaop(sa%l-1)
    bb = shell_dimgaop(sb%l-1)
    bc = shell_dimgaop(sc%l-1)
    do k=1,dc
      hl(1:dleft,1) = h(1:dleft,bc+k)
      call hrr(sa%l,sb%l,q,hl,left)
      do i=1,da
        do j=1,db
          ib(i,j,k) = hl(ba+i,bb+j)!*sa%norm(i)*sb%norm(j)*sc%norm(k)
        end do
      end do
    end do
  end subroutine 

  subroutine integrals_ppi2(ra,rb,sa,sb,ib,nsda)
  ! Evaluate G integrals
  ! Roberto Flores Moreno, 2018
  implicit none
    integer :: nsda
    real(8) :: ra(3),rb(3)
    real(8) :: ib(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb

    integer :: ax,ay,az,ba,bb,bm1,bx,by,bz,da,db,dleft,dright
    integer :: i,ii,im1,im2,j,jj,jm1,k,l,la,lb,ll,lm1,lm2,lt,lx,ly,lz,m,mp1
    real(8) :: f1,f2,f3,f(0:2*nmaxli),zp,zp2,zq,zq2,zw,zw2
    real(8) :: ew,wp,wq,zpw,zqw
    real(8) :: p(3),q(3),w(3)

    integer :: allocs
    real(8),allocatable :: v(:,:,:)

    la = sa%l 
    lb = sb%l
    lt = la + lb 
    dleft = psdaop(la)
    dright = psdaop(lb)
    da = (sa%l+1)*(sa%l+2)/2
    db = (sb%l+1)*(sb%l+2)/2

    allocate(v(dleft,dright,0:lt),stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi2: allocation failed')

    zp = sa%z(1) 
    zp2 = 2.0*zp
    p = ra
    f1 = 2.0*sa%c(1)*((sa%z(1))**nsda)
    zq = sb%z(1) 
    zq2 = 2.0*zq
    q = rb
    zw = zp + zq
    zw2 = 2*zw
    rho = zp*zq/zw
    ew = rho*sum((p-q)**2)
    w = (p*zp + q*zq)/zw
    f2 = f1*sb%c(1)*sqrt(rho/pi)*((pi**2/(zp*zq))**1.5)
    call gammaf(lt,ew,f)
    f(0:lt) = f2*f(0:lt)
    ! s-s-s-s type
    v(1,1,0:lt) = f(0:lt)
    ! vrr right side
    zpw = zp/zw
    do ll=1,lb
      lm1 = ll - 1
      lm2 = max(0,ll-2)
      m = lt - ll
      mp1 = m + 1
      wq = w(3) - q(3)
      ii = gaop(0,0,ll)
      im1 = gaop(0,0,lm1)
      im2 = gaop(0,0,lm2)
      if (lm1.eq.0) then
        v(1,ii,0:m) = wq*v(1,im1,1:mp1) 
      else
        v(1,ii,0:m) = wq*v(1,im1,1:mp1) +       &
     &  (float(lm1)/zq2)*(v(1,im2,0:m)-zpw*v(1,im2,1:mp1))
      end if
      wq = w(2) - q(2)
      do ly=1,ll
        lm1 = ly - 1
        lm2 = max(0,ly-2)
        lz = ll - ly
        ii = gaop(0,ly,lz)
        im1 = gaop(0,lm1,lz)
        im2 = gaop(0,lm2,lz)
        if (lm1.eq.0) then
          v(1,ii,0:m) = wq*v(1,im1,1:mp1) 
        else
          v(1,ii,0:m) = wq*v(1,im1,1:mp1) +     &
     &    (float(lm1)/zq2)*(v(1,im2,0:m)-zpw*v(1,im2,1:mp1))
        end if
      end do
      wq = w(1) - q(1)
      do lx=1,ll
        lm1 = lx - 1
        lm2 = max(0,lx-2)
        do ly=0,ll-lx
          lz = ll - lx - ly
          ii = gaop(lx,ly,lz)
          im1 = gaop(lm1,ly,lz)
          im2 = gaop(lm2,ly,lz)
          if (lm1.eq.0) then
            v(1,ii,0:m) = wq*v(1,im1,1:mp1) 
          else
            v(1,ii,0:m) = wq*v(1,im1,1:mp1) +   &
     &      (float(lm1)/zq2)*(v(1,im2,0:m)-zpw*v(1,im2,1:mp1))
          end if
        end do
      end do
    end do
    ! vrr left side
    zqw = zq/zw
    do ll=1,la
      lm1 = ll - 1
      lm2 = max(0,ll-2)
      m = la - ll
      mp1 = m + 1
      wp = w(3)-p(3)
      ii = gaop(0,0,ll)
      im1 = gaop(0,0,lm1)
      im2 = gaop(0,0,lm2)
      if (lm1.eq.0) then
        v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1) 
      else
        v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*( &
     &        v(im2,1:dright,0:m)-zqw*v(im2,1:dright,1:mp1))
      end if
      do bz=1,lb
        f3 = float(bz)/zw2
        bm1 = bz - 1
        do bx=0,lb-bz
          do by=0,lb-bx-bz
            jj = gaop(bx,by,bz)
            jm1 = gaop(bx,by,bm1)
            v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
          end do
        end do
      end do
      wp = w(2)-p(2)
      do ly=1,ll
        lm1 = ly - 1
        lm2 = max(0,ly-2)
        lz = ll - ly
        ii = gaop(0,ly,lz)
        im1 = gaop(0,lm1,lz)
        im2 = gaop(0,lm2,lz)
        if (lm1.eq.0) then
          v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1)
        else
          v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*( &
     &          v(im2,1:dright,0:m)-zqw*v(im2,1:dright,1:mp1))
        end if
        do by=1,lb
          f3 = float(by)/zw2
          bm1 = by - 1
          do bx=0,lb-by
            do bz=0,lb-bx-by
              jj = gaop(bx,by,bz)
              jm1 = gaop(bx,bm1,bz)
              v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
            end do
          end do
        end do
      end do
      wp = w(1)-p(1)
      do lx=1,ll
        lm1 = lx - 1
        lm2 = max(0,lx-2)
        do ly=0,ll-lx
          lz = ll - lx - ly
          ii = gaop(lx,ly,lz)
          im1 = gaop(lm1,ly,lz)
          im2 = gaop(lm2,ly,lz)
          if (lm1.eq.0) then
            v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1) 
          else
            v(ii,1:dright,0:m) = wp*v(im1,1:dright,1:mp1) + (float(lm1)/zp2)*(&
     &            v(im2,1:dright,0:m)-zqw*v(im2,1:dright,1:mp1))
          end if
          do bx=1,lb
            f3 = float(bx)/zw2
            bm1 = bx - 1
            do by=0,lb-bx
              do bz=0,lb-bx-by
                jj = gaop(bx,by,bz)
                jm1 = gaop(bm1,by,bz)
                v(ii,jj,0:m) = v(ii,jj,0:m)+ f3*v(im1,jm1,1:mp1)
              end do
            end do
          end do
        end do
      end do
    end do
    ba = shell_dimgaop(sa%l-1)
    bb = shell_dimgaop(sb%l-1)
    do i=1,da
      do j=1,db
        ib(i,j) = v(ba+i,bb+j,0)*sa%norm(i)*sb%norm(j)
      end do
    end do
    deallocate(v,stat=allocs)
    if (allocs.gt.0) call file_error('integrals_ppi2: deallocation failed')
  end subroutine 

  subroutine integrals_dipole(ra,rb,sa,sb,ib,nsd)
  ! Evaluate dipole integrals
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    integer :: nsd
    real(8) :: ra(3),rb(3),ib(nmaxnco,nmaxnco,3)
    type(nshell) :: sa,sb

    integer :: ap,ax,ay,az,nb
    type(nshell) :: wsa
    real(8) :: ob(nmaxnco,nmaxnco)

    nb = (sb%l+1)*(sb%l+2)/2

    wsa = sa
    wsa%norm(:) = 1.0
    wsa%l = sa%l + 1
    call integrals_overlap(ra,rb,wsa,sb,ob,nmaxnco,nsd)
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        ap = raop(ax,ay,az)
        ib(ap,1:nb,1) = ob(raop(ax+1,ay,az),1:nb) 
        ib(ap,1:nb,2) = ob(raop(ax,ay+1,az),1:nb) 
        ib(ap,1:nb,3) = ob(raop(ax,ay,az+1),1:nb) 
      end do
    end do

    wsa%l = sa%l
    call integrals_overlap(ra,rb,wsa,sb,ob,nmaxnco,nsd)
    do ax=0,sa%l
      do ay=0,sa%l-ax
        az = sa%l - ax - ay
        ap = raop(ax,ay,az)
        ib(ap,1:nb,1) = ib(ap,1:nb,1) + ra(1)*ob(ap,1:nb) 
        ib(ap,1:nb,2) = ib(ap,1:nb,2) + ra(2)*ob(ap,1:nb) 
        ib(ap,1:nb,3) = ib(ap,1:nb,3) + ra(3)*ob(ap,1:nb) 
        ib(ap,1:nb,:) = ib(ap,1:nb,:)*sa%norm(ap)
      end do
    end do
  end subroutine 

end module

