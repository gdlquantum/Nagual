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
module nmmatrix
! Manipulation of matrices

  use nmfile
  use nmdiag
  use nmstring

  implicit none
  
    public :: matrix_symmetrize
    public :: matrix_diagonalize
    public :: matrix_print
    public :: matrix_print_special
    public :: matrix_print_labeled
    public :: matrix_transpose
    public :: matrix_svd_power
    public :: matrix_multiply
    public :: matrix_transform
    public :: matrix_determinant

    private

contains

  subroutine matrix_print(tape,m,dr,dc,sr,sc,er,ec)
  ! Print out matrix
  ! Roberto Flores-Moreno, Jul 2008
  implicit none
    integer :: dr,dc,sr,sc,er,ec,tape
    real(8) :: m(dr,dc)
    integer :: i,j,nc,c1,c2

    nc = 5
    c2 = sc - 1
    do while (c2.lt.ec)
      c1 = c2 + 1
      c2 = min(c1+nc-1,ec)
      write (tape,'(t2,5i15)') (j,j=c1,c2)
      do i=sr,er
        write (tape,'(t2,5f15.5)') (m(i,j),j=c1,c2)
      end do
    end do
    call file_flush(tape)
  end subroutine 

  subroutine matrix_print_special(tape,m,dr,dc,sr,sc,er,ec)
  ! Print out matrix
  ! Roberto Flores-Moreno, Jul 2008
  implicit none
    integer :: dr,dc,sr,sc,er,ec,tape
    real(8) :: m(dr,dc)
    integer :: i,j,nc,c1,c2

    nc = 13
    c2 = sc - 1
    do while (c2.lt.ec)
      c1 = c2 + 1
      c2 = min(c1+nc-1,ec)
      write (tape,'(t2,13i6)') (j,j=c1,c2)
      do i=sr,er
        write (tape,'(t2,13f6.2)') (m(i,j),j=c1,c2)
      end do
    end do
    call file_flush(tape)
  end subroutine 

  subroutine matrix_print_labeled(tape,m,dr,dc,sr,sc,er,ec,lr,lc)
  ! Print out matrix
  ! Roberto Flores-Moreno, Jul 2008
  implicit none
    character*(*) :: lr(*),lc(*)
    integer :: dr,dc,sr,sc,er,ec,tape
    real(8) :: m(dr,dc)
    integer :: i,j,nc,c1,c2

    nc = 4
    c2 = sc - 1
    do while (c2.lt.ec)
      c1 = c2 + 1
      c2 = min(c1+nc-1,ec)
      write (tape,'(t12,5a15)') (lc(j),j=c1,c2)
      do i=sr,er
        write (tape,'(t2,a5,x,5f15.5)') lr(i),(m(i,j),j=c1,c2)
      end do
    end do
    call file_flush(tape)
  end subroutine 

  subroutine matrix_diagonalize(m,eigv,dm,n)
  ! Diagonalization driver
  ! Roberto Flores Moreno, Jul 2008
  implicit none
    integer :: dm,n
    real(8) :: m(dm,dm),eigv(*)

    call matrix_symmetrize(m,dm,n,'uplow')
    call diag_lapack_dsyev(m,eigv,dm,n)
  end subroutine 

  !>
  !! \f[{\bf C} = {\bf A}\cdot {\bf B}\f]
  !<
  subroutine matrix_multiply(a,b,c,dra,drb,dcb,opa,opb)
  ! Multiply square matrices: c=a.b using BLAS library
  ! Roberto Flores Moreno, Sep 2008
  implicit none
    character :: opa,opb
    integer :: dra,drb,dcb
    real(8) :: a(dra,*),b(drb,*),c(dra,*)

    call dgemm(string_to_uppercase(opa),string_to_uppercase(opb),&
               dra,dcb,drb,1.0,a,dra,b,drb,0.0,c,dra)
  end subroutine

  subroutine matrix_symmetrize(matrix,dmat,n,flag)
  ! Symmetrize matrix
  ! Roberto Flores-Moreno, Dec 2008
  implicit none
    character*(*) :: flag
    integer :: dmat,n
    real(8) :: matrix(dmat,dmat)

    integer :: i,j

    if (flag.eq.'uplow') then
      do i=1,n-1
        do j=i+1,n
          matrix(j,i) = matrix(i,j)
        end do
      end do
    else if (flag.eq.'lowup') then
      do i=1,n-1
        do j=i+1,n
          matrix(i,j) = matrix(j,i)
        end do
      end do
    else 
      call file_error('matrix_symmetrize: unknown flag')
    end if
  end subroutine 

  subroutine matrix_transpose(m,dm)
  ! Transpose a square matrix
  ! Roberto Flores-Moreno, Dec 2008
  implicit none
    integer :: dm
    real(8) :: m(dm,dm)

    integer :: i,j
    real(8) :: hold

    do i=1,dm
      do j=i+1,dm
        hold = m(i,j)
        m(i,j) = m(j,i)
        m(j,i) = hold
      end do
    end do
  end subroutine

  subroutine matrix_transform(m,t,w,dm)
  ! Similarity transformation of m by t: t'mt
  ! Roberto Flores-Moreno, Dec 2008
  implicit none
    integer :: dm
    real(8) :: m(dm,dm),t(dm,dm),w(dm,dm)

    call matrix_multiply(t,m,w,dm,dm,dm,'t','n')
    call matrix_multiply(w,t,m,dm,dm,dm,'n','n')
  end subroutine

  subroutine matrix_svd_power(m,eig,dm,x,tol)
  ! Get negative power of a matrix using singular value decomposition
  ! Roberto Flores-Moreno, Aug 2008
  implicit none
    integer :: dm
    real(8) :: tol,x
    real(8) :: m(dm,dm),eig(*)

    integer :: i,j

    call matrix_diagonalize(m,eig,dm,dm)
    do i=1,dm
      if (eig(i).le.tol) then
        eig(i) = 0.0
        do j=1,dm
          m(j,i) = 0.0
        end do
      else
        do j=1,dm
          if (x.lt.0) then
            m(j,i) = m(j,i)/sqrt(eig(i))**abs(x)
          else
            m(j,i) = m(j,i)*sqrt(eig(i))**abs(x)
          end if
        end do
      end if
    end do
    m = matmul(m,transpose(m))
  end subroutine

  subroutine matrix_invert_by_diag(m,eig,dm)
  ! Get negative power of a matrix using singular value decomposition
  ! Roberto Flores-Moreno, Aug 2008
  implicit none
    integer :: dm
    real(8) :: m(dm,dm),eig(*)

    integer :: i,j

    call matrix_diagonalize(m,eig,dm,dm)
    do i=1,dm
      do j=1,dm
        m(j,i) = m(j,i)/eig(i)**0.5
      end do
    end do
    m = matmul(m,transpose(m))
  end subroutine

  real(8) recursive function matrix_determinant(m,dm) result( det )
  ! Evaluates matrix determinant
  ! Roberto Flores-Moreno (Oct 2014)
  implicit none
    integer :: dm
    real(8) :: m(dm,dm)

    integer :: i,j,jj,k,kk
    real(8) :: detm,signus

    integer :: allocs,dminor
    real(8),allocatable :: minor(:,:)

    if (dm.eq.2) then
      det = m(1,1)*m(2,2) - m(1,2)*m(2,1)
    else if (dm.gt.2) then
      dminor = dm - 1
      allocate(minor(dminor,dminor),stat=allocs)
      if (allocs.gt.0) call file_error('matrix_determinant: allocation failed')

      signus = -1.0
      det = 0.0
      do i=1,dm
        signus = -signus
        do j=2,dm
          jj = j - 1
          do k=1,dm
            if (k.ne.i) then
              kk = k
              if (k.gt.i) kk = kk - 1
              minor(jj,kk) = m(j,k)
            end if
          end do
        end do
        detm = matrix_determinant(minor,dminor)
        det = det + m(1,i)*signus*detm
      end do

      deallocate(minor,stat=allocs)
      if (allocs.gt.0) call file_error('matrix_determinant: deallocation failed')
    else if (dm.eq.1) then
      det = m(1,1)
    else if (dm.lt.1) then
      det = 0.0
    end if
  end function

end module
