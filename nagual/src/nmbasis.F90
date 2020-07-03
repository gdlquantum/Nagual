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
module nmbasis
  
  use nmparameter
  use nmtypes
  use nmfile
  use nmset
  use nmshell
  use nmmolecule
  use nmmath

  implicit none

    public :: basis_initialize
    public :: basis_nfun
    public :: basis_restricted_nfun
    public :: basis_nfun_spherical
    public :: basis_copy
    public :: basis_print
    public :: basis_lmax
    public :: basis_radial_domains

    real(8), public :: ctostm(2*nmaxlb+1,(nmaxlb+1)*(nmaxlb+2)/2,0:nmaxlb)
    real(8), public :: ptostm(nmaxnco,nmaxnco*nmaxli,0:nmaxli)

    private

contains

  subroutine basis_initialize(b)
  ! Initialize basis
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nbasis) :: b

    call basis_normalize(b)
    call basis_radial_domains(b)
    call basis_transformation_factors
  end subroutine

  subroutine basis_normalize(b)
  ! Normalization
  ! Roberto Flores-Moreno, 2018
  implicit none 
    type(nbasis) :: b

    integer :: i,iset,ishell,j,lx,ly,lz
    real(8) :: factor,norm
    type(nshell) :: s

    do iset=1,b%nsets
      do ishell=1,b%set(iset)%nshell
        s = b%set(iset)%shell(ishell)
        norm = 2.0**s%l*(2.0/pi)**0.75
        do i=1,s%k
          s%c(i) = norm*(s%z(i)**(2*s%l+3))**0.25*s%c(i)
        end do
        factor = pi**1.5/2**s%l
        norm = 0.0
        do i=1,s%k
          do j=1,s%k
            norm = norm + s%c(i)*s%c(j)*factor/(s%z(i) + s%z(j))**(s%l+1.5)
          end do
        end do
        do lx=0,s%l
          do ly=0,s%l-lx
            lz = s%l - lx - ly
            s%norm(raop(lx,ly,lz)) = math_factoriald(2*lx-1)*   &
     &              math_factoriald(2*ly-1)*math_factoriald(2*lz-1)
            s%norm(raop(lx,ly,lz)) = 1.0/sqrt(norm*s%norm(raop(lx,ly,lz)))
          end do
        end do
        b%set(iset)%shell(ishell) = s
      end do
    end do
  end subroutine

  integer function basis_nfun(b)
  ! Get number of basis functions 
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    type(nbasis) :: b

    integer :: iset,ishell,l

    basis_nfun = 0
    do iset=1,b%nsets
      do ishell=1,b%set(iset)%nshell
        l = b%set(iset)%shell(ishell)%l
        basis_nfun = basis_nfun + (l+1)*(l+2)/2
      end do
    end do
  end function 

  integer function basis_restricted_nfun(b,ptr)
  ! Get number of basis functions on given set of atoms
  ! Roberto Flores-Moreno, 2018
  implicit none
    logical :: ptr(*)
    type(nbasis) :: b

    integer :: iset,ishell,l

    basis_restricted_nfun = 0
    do iset=1,b%nsets
      if (ptr(iset)) then
        do ishell=1,b%set(iset)%nshell
          l = b%set(iset)%shell(ishell)%l
          basis_restricted_nfun = basis_restricted_nfun + (l+1)*(l+2)/2
        end do
      end if
    end do
  end function 

  integer function basis_nfun_spherical(b)
  ! Get number of spherical  basis functions
  ! Roberto Flores-Moreno, 2013
  implicit none
    type(nbasis) :: b

    integer :: iset,ishell,l

    basis_nfun_spherical = 0
    do iset=1,b%nsets
      do ishell=1,b%set(iset)%nshell
        l = b%set(iset)%shell(ishell)%l
        basis_nfun_spherical = basis_nfun_spherical + 2*l + 1
      end do
    end do
  end function 

  subroutine basis_copy(src,dst)
  ! Copy basis
  ! Roberto Flores Moreno, 2018
  implicit none
    type(nbasis) :: src,dst

    integer :: iset

    dst%nsets = src%nsets

    do iset=1,src%nsets
      call set_copy(src%set(iset),dst%set(iset))
    end do
  end subroutine 

  subroutine basis_print(tape,b)
  ! Print molecule basis to a file
  ! Roberto Flores Moreno, 2018
  implicit none
    integer :: tape
    type(nbasis) :: b

    integer :: iset

    if (b%nsets.le.0) return

    do iset=1,b%nsets
      write(tape,'(/,"%%%%%%%%%%%%%",/,t2,i6)') iset
      call set_print(tape,b%set(iset))
    end do
    write(tape,'(/,"%%%%%%%%%%%%%")') 
    call file_flush(tape)
  end subroutine 

  integer function basis_lmax(b)
  ! Get maximum l value
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nbasis) :: b

    integer :: iset,l

    basis_lmax = 0
    do iset=1,b%nsets
      l = set_lmax(b%set(iset))
      if (l.gt.basis_lmax) basis_lmax = l
    end do
  end function

  subroutine basis_transformation_factors
  ! Basis/ecp transformation factors
  ! H. B. Schlegel and M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995).
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    integer :: expo,i,ica,ic1,ic2,is,isa,j,k,l,l1,l2,lx,lx1,lx2
    integer :: ly,ly1,ly2,lz,lz1,lz2,m,ma,nc
    real(8) :: pi,s,s1,s2,xyz
    real(8) :: fac(0:2*nmaxlb)

    ! Factorials
    fac(0) = 1.0
    do i=1,2*nmaxlb
      fac(i) = fac(i-1)*float(i)
    end do

    do l=0,nmaxlb
      ica = 0
      do lx=l,0,-1
        do ly=l-lx,0,-1
          lz = l - lx - ly
          ica = ica + 1
          do m=-l,l
            isa = l + m + 1
            ma = abs(m)
            j = lx + ly - ma
            if ((j.ge.0).and.(mod(j,2).eq.0)) then
              j = j/2
              s1 = 0.0
              do i=0,(l-ma)/2
                s2 = 0.0
                do k=0,j
                  if (((m.lt.0).and.(mod(abs(ma-lx),2).eq.1)).or.       &
     &                ((m.gt.0).and.(mod(abs(ma-lx),2).eq.0))) then
                    expo = (ma - lx + 2*k)/2
                    s = (-1.0)**expo*sqrt(2.0)
                  else if ((m.eq.0).and.(mod(lx,2).eq.0)) then
                    expo = k - lx/2
                    s = (-1.0)**expo
                  else
                    s = 0.0
                  end if
                  if ((lx-2*k.ge.0).and.(lx-2*k.le.ma)) then
                    s2 = s2 + (fac(j)/(fac(k)*fac(j-k)))*               &
     &              (fac(ma)/(fac(ma-(lx-2*k))*fac(lx-2*k)))*s
                  end if
                end do
                if (j.le.i) then
                  s1 = s1 + (fac(l)/(fac(l-i)*fac(i)))*                 &
     &               (fac(i)/(fac(i-j)*fac(j)))*(-1.0)**i*              &
     &               fac(2*l-2*i)/fac(l-ma-2*i)*s2
                end if
              end do
              ctostm(isa,ica,l) = sqrt((fac(2*lx)*fac(2*ly)*fac(2*lz)*  &
     &        fac(l)*fac(l-ma))/(fac(lx)*fac(ly)*fac(lz)*               &
     &        fac(2*l)*fac(l+ma)))*s1/(2.0**l*fac(l))
            else
              ctostm(isa,ica,l) = 0.0
            end if
          end do
        end do
      end do
    end do

    ! Initialize PTOSTM 
    do l=0,nmaxli
      nc = ((l+1)*(l+2))/2
      do is=1,2*l+1
        do ic1=1,nc
          ptostm(is,ic1,l) = 0.0
        end do
      end do
    end do

    pi = 2.0*acos(0.0)

    ! Build the transformation matrix for the transformation
    ! of a polynomial in a sum of spherical harmonics (PTOSTM) 
    do l1=0,nmaxli
      do lx1=0,l1
        do ly1=0,l1-lx1
          lz1 = l1 - lx1 - ly1
          ic1 = gaop(lx1,ly1,lz1)
          do l2=0,l1
            s1 = 4.0*pi*math_factoriald(2*l2+1)
            do is=1,2*l2+1
              xyz = 0.0
              do lx2=0,l2
                do ly2=0,l2-lx2
                  lz2 = l2 - lx2 - ly2
                  ic2 = raop(lx2,ly2,lz2)
                  lx = lx1 + lx2
                  ly = ly1 + ly2
                  lz = lz1 + lz2
                  if ((mod(lx,2).eq.0).and.(mod(ly,2).eq.0).and.&
                      (mod(lz,2).eq.0)) then
                    s2 = math_factoriald(lx-1)*math_factoriald(ly-1)*&
                         math_factoriald(lz-1)/math_factoriald(l1+l2+1)
                    s = s1/(math_factoriald(2*lx2-1)*math_factoriald(2*ly2-1)*&
                            math_factoriald(2*lz2-1))
                    xyz =  xyz + ctostm(is,ic2,l2)*sqrt(s)*s2
                  end if
                end do
              end do
              ptostm(is,ic1,l2) = xyz
            end do
          end do
        end do
      end do
    end do
  end subroutine

  subroutine basis_radial_domains(b)
  ! Evaluate shell radius
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nbasis) :: b

    integer :: ba,i,iset,ishell,na
    real(8) :: alpha,d,r
    type(nshell) :: sa

    ba = 0
    do iset=1,b%nsets
      b%set(iset)%radius = 0.0
      b%set(iset)%ll = ba + 1
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = shell_nco(sa%l)
        alpha = sa%z(1)
        d = abs(sa%c(1))
        do i=2,sa%k
          if (sa%z(i).lt.alpha) then
            alpha = sa%z(i)
            d = abs(sa%c(i))
          end if
        end do
        r = shell_gto_radius(alpha,d,sa%l,1.0e-10)
        b%set(iset)%shell(ishell)%radius = r
        b%set(iset)%radius = max(r,b%set(iset)%radius)
        ba = ba + na
      end do
      b%set(iset)%ul = ba
    end do
  end subroutine

end module
