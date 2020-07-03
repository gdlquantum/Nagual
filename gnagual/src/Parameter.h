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
#ifndef GN_PARAMETER_H
#define GN_PARAMETER_H
 
// Tolerances
#define GN_TOL_NUM  1.0e-14

// Strings
#define GN_MAX_STR_SIZE 256

// Orbital calculations
#define GN_MAX_ATOMIC_NUMBER 103+2
#define GN_MAX_CON      20
#define GN_MAX_L_DER     0
#define GN_MAX_L_BASIS   4
#define GN_MAX_L_ECPS    4
#define GN_MAX_L       GN_MAX_L_BASIS + GN_MAX_L_DER
#define GN_MAX_L_I     GN_MAX_L + GN_MAX_L_BASIS
#define GN_MAX_NCO     ((GN_MAX_L+1)*(GN_MAX_L+2))/2
#define GN_MAX_TAB_GAM 120
#define GN_MAX_NRQP    200
#define GN_MAX_NSPEC 1
#define GN_MAX_NSPIN 2

// Symmetry
#define GN_SYMMETRY_TOL 1.0e-5

// Grid
#define GN_GS  5.0e-4       //  5.0 for eem, 5.0e-4 for nuclear structure

#endif // GN_PARAMETER_H
