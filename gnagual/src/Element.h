/*
! This file is part of GNagual software.
!
!    GNagual is free software: you can redistribute it and/or modify
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
! R. Flores-Moreno, 
! GNagual 0.7, Guadalajara Jal., Mexico (2019)
!
!###################################################
!   GNagual: Graphics for Nagual
!   Copyright (C) 2003-2019 GNagual developers.
!
!   List of authors 
!
!   R. Flores-Moreno            roberto.floresmoreno.qt@gmail.com
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
#ifndef GN_ELEMENT_H
#define GN_ELEMENT_H

#include <Parameter.h>

extern const char*  ELEMENT_SYMBOL[GN_MAX_ATOMIC_NUMBER+1];
extern const double ELEMENT_COV_R[GN_MAX_ATOMIC_NUMBER+1];
extern const double ELEMENT_VDW_R[GN_MAX_ATOMIC_NUMBER+1];
extern const double ELEMENT_MASS[GN_MAX_ATOMIC_NUMBER+1];
extern const char*  ELEMENT_NAME[GN_MAX_ATOMIC_NUMBER+1];
extern double ELEMENT_COLOR[GN_MAX_ATOMIC_NUMBER+1][3];

#endif
