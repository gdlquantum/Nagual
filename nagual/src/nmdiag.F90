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
module nmdiag
! Matrix diagonalizations in different ways

  use nmfile

  implicit none

  public :: diag_lapack_dsyev

  private

contains

  subroutine diag_lapack_dsyev(m,eigv,dm,n)
  ! Diagonalize a matrix completely using LAPACK dsyev
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: dm,n
    real(8) :: m(dm,dm),eigv(*)

    integer :: allocs,drw1,ierr
    real(8),allocatable :: rw1(:)
    
    drw1 = max(1,3*dm-1)
    allocate(rw1(drw1),stat=allocs)
    if (allocs.gt.0) call file_error('diag_lapack_dsyev: allocation')

    call dsyev('V','U',n,m,dm,eigv,rw1,drw1,ierr)
    if (ierr.ne.0) then
      call file_error('diagonalization failed DSYEV')
    end if

    deallocate(rw1,stat=allocs)
    if (allocs.gt.0) call file_error('diag_lapack_dsyev: deallocation')
  end subroutine 

end module
