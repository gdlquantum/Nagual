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
module nmtime
! Time control utilities

  use nmtypes
  use nmstring
  use nmfile

  implicit none

    public :: start_timer
    public :: stop_timer
    public :: continue_timer
    public :: print_timer

    private

contains

  subroutine start_timer(timer,title)
  ! Initialize timer
  ! Roberto Flores-Moreno, Sep 2008
  implicit none
    type(ntimer) :: timer
    character*(*) :: title
   
    timer%title = title
    timer%val = 0.0
    call continue_timer(timer)
  end subroutine

  subroutine stop_timer(timer)
  ! Stop timer
  ! Roberto Flores-Moreno, Sep 2008
  implicit none
    type(ntimer) :: timer
   
    real(4) :: t

    timer%stat = 0
    t = gettime()
    timer%val = timer%val + t - timer%when
    timer%when = t
  end subroutine

  subroutine continue_timer(timer)
  ! Continue timer
  ! Roberto Flores-Moreno, Sep 2008
  implicit none
    type(ntimer) :: timer
   
    if (timer%stat.eq.1) return
    timer%stat = 1
    timer%when = gettime()
  end subroutine

  subroutine print_timer(timer,tape)
  ! Print timer value
  ! Roberto Flores-Moreno, Sep 2008
  implicit none
    type(ntimer) :: timer
    integer :: tape

    real(8) :: nh,nm,ns
    real(4) t
   
    if (timer%stat.eq.1) then
      t = gettime()
      timer%val = timer%val + t - timer%when
      timer%when = t
    end if
    nh = timer%val/3600
    nm = timer%val/60
    ns = timer%val

    if (nh.gt.1) then
      write(tape,'(/,t2,a50,x,f10.2,x,"hours")') timer%title,nh
    else if (nm.gt.1) then
      write(tape,'(/,t2,a50,x,f10.2,x,"minutes")') timer%title,nm
    else
      write(tape,'(/,t2,a50,x,f10.0,x,"seconds")') timer%title,timer%val
    end if

    call file_flush(tape)
  end subroutine

  real(4) function gettime()
  ! This function takes time as real value in seconds
  ! Roberto Flores-Moreno, Sep 2008
  implicit none
    real(4) :: etime,tv(2)

    tv(1:2) = 0.0
    gettime = etime(tv)
  end function

end module
