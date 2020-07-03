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
module nmgrid
! Grid

  use nmtypes
  use nmfile
  use nmbecke
  use nmlebedev
  use nmbasis
  use nmset
  use nmshell

  implicit none

    public :: grid_initialize
    public :: grid_basis

    integer,public :: ngpc(nmaxatom)

    private 

contains

  subroutine grid_initialize(sys)
  ! Initialize grid
  ! Roberto Flores-Moreno, 2019
  implicit none
    type(nsystem) :: sys

    integer :: is

    do is=1,sys%ns
      call basis_radial_domains(sys%auxis(is))
    end do

    call grid_fixed(sys,1)
  end subroutine

  subroutine grid_fixed(sys,is)
  ! Fixed grid
  ! Roberto Flores-Moreno, 2020
  implicit none
    integer :: is
    type(nsystem) :: sys

    integer :: i,iatom,irp,iset
    real(8) :: r0(3)

    integer :: allocs,dap,drp
    real(8),allocatable :: rp(:),rw(:)
    real(8),allocatable :: xa(:),ya(:),za(:),wa(:)
    real(8),allocatable :: g(:,:),gbuf(:,:)

    print *, "Generating grid "
    dap = 302
    drp = 75
    print *, "Number of radial points: ",drp
    print *, "Number of angular points: ",dap
    if (dap*drp.gt.nmaxgpa) call file_error('grid_fixed too many points')
    allocate(rp(drp),rw(drp),g(4,dap),gbuf(4,nmaxgpa),xa(dap),&
             ya(dap),za(dap),wa(dap),stat=allocs)
    if (allocs.gt.0) call file_error('grid_fixed, allocation failed')

    call becke_radial_quadrature(rp,rw,drp)
    rw(1:drp) = 8.0*acos(0.0)*rw(1:drp)*rp(1:drp)**2

    call lebedev(xa,ya,za,wa,dap)

    do iset=1,sys%auxis(is)%nsets
      ngpc(iset) = 0
      iatom = iset
      r0(1:3) = sys%mol%atom(iatom)%pos(1:3)
      do irp=1,drp
        do i=1,dap
          g(1,i) = xa(i)*rp(irp) + r0(1)
          g(2,i) = ya(i)*rp(irp) + r0(2)
          g(3,i) = za(i)*rp(irp) + r0(3)
          g(4,i) = wa(i)*rw(irp)*becke_weight(sys%mol,iset,g(1,i))
          gbuf(1:4,ngpc(iset)+i) = g(1:4,i)
        end do
        ngpc(iset) = ngpc(iset) + dap
      end do
      call file_grid(gbuf,4*nmaxgpa,iset,'write')
    end do

    deallocate(rp,rw,g,gbuf,xa,ya,za,wa,stat=allocs)
    if (allocs.gt.0) call file_error('grid_fixed, deallocation failed')
  end subroutine

  subroutine grid_vicinity(m,b,ref,n,nn,ds)
  ! Build system with a given atom and its neighbors
  ! Roberto Flores-Moreno, 2019
  implicit none
    integer :: ds,nn,ref
    integer :: n(*)
    type(nmolecule) :: m
    type(nbasis) :: b

    integer :: iatom,iset,ishell,l
    real(8) :: d0,d,dcut
    real(8) :: r0(3),ra(3),rr(3)

    iatom = ref
    r0(1:3) = m%atom(iatom)%pos(1:3)
    d0 = sqrt(b%set(ref)%radius)

    ds = 0
    nn = 0
    do iset=1,b%nsets
      iatom = iset
      ra(1:3) = m%atom(iatom)%pos(1:3)
      rr(1:3) = ra(1:3) - r0(1:3)
      d = sqrt(sum(rr(1:3)*rr(1:3)))
      dcut = sqrt(b%set(iset)%radius) + d0
      if (d.le.dcut) then
        nn = nn + 1
        n(nn) = iset
        do ishell=1,b%set(iset)%nshell
          l = b%set(iset)%shell(ishell)%l
          ds = ds + shell_nco(l)
        end do
      end if
    end do
  end subroutine

  subroutine grid_basis(m,basis,g,bas,dbas,ngp)
  ! Evaluates auxiliary functions
  ! Roberto Flores-Moreno, 2009, 2018
  implicit none
    integer :: dbas,ngp
    real(8) :: g(4,*),bas(dbas,*)
    type(nmolecule) :: m
    type(nbasis) :: basis

    integer :: ba,i,iatom,igp,iset,irp,ishell,l,lmax,lx,ly,lz,na
    real(8) :: expf,r2,z,zp
    real(8) :: ra(3),p(3,0:nmaxl)
    type(nshell) :: sa

    bas(1:dbas,1:ngp) = 0.0

    p(:3,0) = 1.0
    do igp=1,ngp
      ba = 0
      do iset=1,basis%nsets
        iatom = iset
        ra(1:3) = g(1:3,igp) - m%atom(iatom)%pos(1:3)
        r2 = sum(ra*ra)
        if (r2.lt.basis%set(iset)%radius) then
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
                bas(ba+i,igp) = expf*sa%norm(i)*p(1,lx)*p(2,ly)*p(3,lz)
              end do
            end do
            ba = ba + na
          end do
        else
          do ishell=1,basis%set(iset)%nshell
            sa = basis%set(iset)%shell(ishell)
            na = shell_nco(sa%l)
            ba = ba + na
          end do
        end if
      end do
    end do
  end subroutine

end module

