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
module nmsort

  implicit none

    public :: sort_vector

    private

contains

  subroutine sort_vector(r,ptr,n)
  ! Sort vector 
  ! W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery,
  ! Numerical Recipes, 2nd Edition, (Cambridge University Press)
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: n
    integer :: ptr(*)
    real(8) :: r(*)

    integer :: i,indx,j,sel,tmpptr
    real(8) :: tmp

    if (n.lt.2) return
    do i=1,n
      ptr(i) = i
    end do

    indx = n/2+1
    sel= n
 10 continue
      if (indx.gt.1) then
        indx = indx-1
        tmp = r(indx)
        tmpptr = ptr(indx)
      else
        tmp = r(sel)
        tmpptr = ptr(sel)
        r(sel) = r(1)
        ptr(sel) = ptr(1)
        sel = sel - 1
        if (sel.eq.1) then
          r(1) = tmp
          ptr(1) = tmpptr
          go to 1000
        end if
      end if
      i = indx
      j = 2*indx
 20   continue
      if (j.le.sel) then
        if (j.lt.sel) then
          if (r(j).lt.r(j+1)) j = j + 1
        end if
        if (tmp.lt.r(j)) then
          r(i) = r(j)
          ptr(i) = ptr(j)
          i = j
          j = 2*j
        else
          j = sel + 1
        end if
      go to 20
      end if
      r(i) = tmp
      ptr(i) = tmpptr
    go to 10
1000 continue
  end subroutine

end module
