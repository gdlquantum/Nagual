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
#ifndef GN_MATH_H
#define GN_MATH_H

#include <math.h>

#define GN_PI M_PI //2.0*acos(0.0)
#define GN_ABS(a) ((a)>=0.0?(a):-(a))
#define GN_MOD(a,b) (((int)(a))-((int)(b))*(((int)(a))/((int)(b)))) 

#define GN_MAX(a,b) ((a)>(b)?(a):(b))
#define GN_MIN(a,b) ((a)<(b)?(a):(b))

extern double Factorial(int);
extern double DoubleFactorial(int);
extern double Noverk(int,int);

#endif
