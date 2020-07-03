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
module nmcfvwn
! VWN correlation functionals

  use nmparameter
  use nmmath

  implicit none

    public :: cfvwn_energy
    public :: cfvwn_potential
    public :: cfvwn_kernel
    public :: cfvwn_kernel2

    private 
      real(8) :: vwns,vwnf2z,r13,r23,r43,r16,r56,r49
      parameter ( vwns = 2.0*(2.0**(1.0/3.0) - 1.0), &
                  vwnf2z = 4.0/(9.0*(2.0**(1.0/3.0) - 1.0)),&
                  r13 = 1.0/3.0, r23 = 2.0/3.0, r43 = 4.0/3.0,&
                  r16 = 1.0/6.0, r56 = 5.0/6.0, r49 = 4.0/9.0 )

      real(8) :: vwna1,vwnb1,vwnc1,vwnx01
      real(8) :: vwna2,vwnb2,vwnc2,vwnx02
      real(8) :: vwna3,vwnb3,vwnc3,vwnx03
      parameter (vwna1=0.0621814,   vwna2=0.0310907,   vwna3=-0.0337737,  &
                 vwnb1=3.7274400,   vwnb2=7.0604200,   vwnb3=1.1310710,   &
                 vwnc1=12.9352000,  vwnc2=18.0578000,  vwnc3=13.0045000,  &
                 vwnx01=-0.1049800, vwnx02=-0.3250000, vwnx03=-0.0047584)

contains
! S. Vosko, L. Wilk, M. Nusair, 
! Can. J. Phys. 58, 1200 (1980)

  real(8) function padevwn(x,a,b,c,x0)
  implicit none
    real(8) :: a,b,c,x,x0

    real(8) :: den,q,x0s,xs

    q = sqrt(4.0*c-b*b)
    xs = x*x
    x0s = x0*x0
    den = (x - x0)
    den = den*den
    padevwn = a*(log(xs/den) - ((x0s + c)*log((xs+b*x+c)/den) +         &
              2.0*b*(x0s-c)*atan(q/(2.0*x+b))/q)/(x0s+b*x0+c))
  end function

  real(8) function dpadevwn(x,a,b,c,x0)
  implicit none
    real(8) :: a,b,c,x,x0

    real(8) :: xs
  ! DPADE(x) = -(x/6) p'(x) 
    xs = x*x
    dpadevwn = r13*a*((1.0+b/(x-x0))*(xs/(xs+b*x+c))-1.0)
  end function

  real(8) function d2padevwn(x,a,b,c,x0)
  implicit none
    real(8) :: a,b,c,x,x0

    real(8) :: xs,dx,px
  ! D2PADE(x) = (x/6)^2 p''(x)
    xs = x*x
    dx = x-x0
    px = xs+b*x+c
    d2padevwn= 1.0/xs-1.0/px*(b*x/(dx**2)-(1.0+b/dx)*((c-xs)/px))
    d2padevwn = -2.0*a*d2padevwn*xs/36.0
  end function

  real(8) function d3padevwn(x,a,b,c,x0)
  implicit none
    real(8) :: a,b,c,x,x0

    real(8) :: dx,px,px2,xfac,x3

    x3 = x*x*x
    dx = x-x0
    px = x**2+b*x+c
    px2 = px*px
    xfac = -2.0*x3 + 6.0*x*c + 2.0*b*c
    d3padevwn = 2.0/x3 + 1.0/px*(xfac/px2*(1.0 + b/dx) + &
               (1.0 - 2.0*x/dx + (c - x*(3.0*x + b))/px)*b/dx**2)
    d3padevwn = 2.0*a*d3padevwn
    d3padevwn = -x3*d3padevwn/216.0
  end function


  real(8) function srs(rhot)
  implicit none
    real(8) :: rhot

    srs = (0.75/(pi*rhot))**r16
  end function

  subroutine cfvwn_energy(rhoa,rhob,ec)
  ! Evaluates VWN correlation energy functional
  ! Roberto Flores-Moreno, 2010, 2018
  implicit none
    real(8) :: rhoa,rhob,ec

    integer :: i
    real(8) :: rhot,s,ep,ef,ez,sp,zs,zs2,zs4

    rhot = rhoa + rhob
    if (rhot.ne.0.0) then
      s = srs(rhot)
      ep = padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
      zs = (rhoa - rhob)/rhot
      ec = ep
    else
      zs = 0.0
      ec = 0.0
    end if
    if (abs(zs).gt.ntolnum) then
      sp = ((1.0+zs)**r43+(1.0-zs)**r43 - 2.0)/vwns
      ef = padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      ez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)
      zs2 = zs*zs
      zs4 = zs2*zs2
      ec = ec + sp*((ef-ep)*zs4 + ez*(1.0-zs4)/vwnf2z)
    end if
    ec = 0.5*ec*rhot
  end subroutine

  subroutine cfvwn_potential(rhoa,rhob,vca,vcb)
  ! Evaluates VWN correlation potential
  ! Roberto Flores-Moreno, 2010, 2018
  implicit none
    real(8) :: rhoa,rhob,vca,vcb

    integer :: i
    real(8) :: def,defpz,dep,dez,dsp,ef,efpz,ep,ez,rhot
    real(8) :: s,sp,vcfp,vcpol,zs,zs3

    rhot = rhoa + rhob
    if (rhot.gt.0.0) then
      s = srs(rhot)
      ep = padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
      dep = dpadevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
      vca = ep + dep
      vcb = vca
      zs = (rhoa - rhob)/rhot
    else
      vca = 0.0
      vcb = 0.0
      zs = 0.0
    end if
    if (abs(zs).gt.ntolnum) then
      sp = ((1.0+zs)**r43+(1.0-zs)**r43 - 2.0)/vwns
      dsp = r43*((1.0+zs)**r13-(1.0-zs)**r13)/vwns
      zs3 = zs*zs*zs
      ef = padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      ez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
      efpz = zs3*(ef - ep - ez)
      def = dpadevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      dez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
      defpz = zs3*(def - dep - dez)
      vcfp = sp*(zs*(efpz + defpz) + ez + dez)
      vcpol = dsp*(zs*efpz + ez) + 4.0*sp*efpz
      vca = vca + vcfp + (1.0 - zs)*vcpol
      vcb = vcb + vcfp - (1.0 + zs)*vcpol
    end if
    vca = 0.5*vca
    vcb = 0.5*vcb
  end subroutine

  subroutine cfvwn_kernel(rhoa,rhob,faa,fab,fbb)
  ! Evaluates VWN correlation kernel
  ! Roberto Flores-Moreno, 2010, 2018
  implicit none
    real(8) :: rhoa,rhob
    real(8) :: faa,fab,fbb

    integer :: i,j,k
    real(8) :: vcfp,vcpol
    real(8) :: czs,zs,zs3,s,rhot
    real(8) :: ep,ef,efpz,ez,sp
    real(8) :: dep,def,defpz,dez,dsp
    real(8) :: d2ep,d2ef,d2efpz,d2ez,d2sp

    rhot = rhoa + rhob
    if (rhot.ne.0.0) then
      zs = (rhoa - rhob)/rhot
      s = srs(rhot)

      d2sp = 0.0
      if (abs(1.0+zs).gt.ntolnum) then
        d2sp = d2sp + (1.0 + zs)**(-r23)
      end if
      if (abs(1.0-zs).gt.ntolnum) then
        d2sp = d2sp + (1.0 - zs)**(-r23)
      end if
      d2sp = r49*d2sp/vwns

      ez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
      dep = dpadevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
      d2ep = d2padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)

      faa = d2ep + r56*dep + d2sp*ez
      fbb = faa
      fab = d2ep + r56*dep - d2sp*ez
    else
      zs = 0.0
      faa = 0.0
      fab = 0.0
      fbb = 0.0
      return
    end if

    if (abs(zs).gt.ntolnum) then
      sp = ((1.0+zs)**r43+(1.0-zs)**r43 - 2.0)/vwns
      dsp = r43*((1.0+zs)**r13-(1.0-zs)**r13)/vwns
      ep = padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
      ef = padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      def = dpadevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      dez = dpadevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
      d2ef = d2padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
      d2ez = d2padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
      zs3 = zs**3
      czs = 1.0 - zs
      efpz = zs3*(ef - ep - ez)
      defpz = zs3*(def - dep - dez)
      d2efpz = zs3*(d2ef - d2ep - d2ez)
      faa = faa + sp*(12.0*efpz/zs*czs**2 - defpz*(zs*(7.0+0.5*r13)-8.0) + &
            r56*dez+zs*d2efpz+d2ez) + 2.0*dsp*czs*(zs*defpz + dez + 4.0* &
            czs*efpz)+ d2sp*((ez + zs*efpz)*czs**2-ez)
      fbb = faa + 16.0*sp*(3.0*efpz - defpz) + 4.0*dsp*(zs*(8.0*efpz-defpz) - &
            dez) + 4.0*d2sp*zs*(ez+zs*efpz)
      fab = faa - 8.0*sp*(3.0*czs*efpz/zs + defpz) - 2.0*dsp*(8.0*efpz*czs + &
            zs*defpz + dez) - 2.0*d2sp*czs*(ez + zs*efpz)
    end if

    faa = 0.5*faa/rhot
    fbb = 0.5*fbb/rhot
    fab = 0.5*fab/rhot
  end subroutine

  subroutine cfvwn_kernel2(rhoa,rhob,gaaa,gaab,gabb,gbbb)
  ! Evaluates VWN correlation second kernel
  ! Roberto Flores-Moreno, 2010, 2018
  implicit none
    real(8) :: rhoa,rhob,gaaa,gaab,gabb,gbbb

    real(8) :: dep,d2ep,d3ep,rhot,s

    rhot = rhoa + rhob
    s = srs(rhot)
    dep = dpadevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
    d2ep = d2padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
    d3ep = d3padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
    gaaa = -0.5*((35.0/36.0)*dep - 0.5*d2ep + d3ep)/rhot**2
    gaab = gaaa
    gabb = gaaa
    gbbb = gaaa
  end subroutine

end module 
