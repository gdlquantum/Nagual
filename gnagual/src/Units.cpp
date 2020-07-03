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
#include <Units.h>
#include <Math.h>

// === Energy ===

// Conversion factor taken from ?
double HartreeToeV(double a) { return ((a)*27.21138466); };
double eVToHartree(double a) { return ((a)/27.21138466); };

// Conversion factor taken from the book
// "Physical Chemistry: A Molecular Approach", D. A. McQuarrie and J. D. Simon,
// University Science Books (Sausalito CA, 1997).
double HartreeToJMol(double a) { return ((a)*2625500.0); };
double JMolToHartree(double a) { return ((a)/2625500.0); };

// === Length ===

// Conversion factor taken from deMon2k 4.0.8
double AngstromToBohr(double a) { return ((a)*1.8897261349309467); }; 
double BohrToAngstrom(double a) { return ((a)/1.8897261349309467); };

// === Angles ===
double DegreeToRadian(double deg) { return GN_PI*deg/180.0; }; 
double RadianToDegree(double rad) { return 180.0*rad/GN_PI; };

