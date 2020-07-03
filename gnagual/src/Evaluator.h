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
#ifndef GN_EVALUATOR_H
#define GN_EVALUATOR_H

#include <vector>

using namespace std;

class System;
class Spin;
class Set;
class Matrix;

extern void GetGTO(vector<Set*>,int,double*,double*,double*,double*,
            int);
extern void GetBasis(vector<Set*>,double,double,double,double*);
extern void GetOrbital(Spin*,int,double*,double*,double*,double*,int);
extern void GetDensityLike(Matrix*,vector<Set*>,double*,double*,double*,
            double*,int);
extern void GetSpeciesDensity(Spin*,int,double*,double*,double*,double*,int);
extern void GetDensity(System*,int,double*,double*,double*,double*,int);
extern void GetChargeDensity(System*,int,double*,double*,double*,double*,int);

extern void GetShape(System*,int,double*,double*,double*,double*,int);
extern void GetShannon(System*,int,double*,double*,double*,double*,int);
extern void GetElectronicSpinDensity(System*,int,double*,double*,double*,
            double*,int);

extern void GetGTODerivatives(vector<Set*>,int,double*,double*,double*,
            double*,double*,double*,int);
extern void GetBasisDerivatives(vector<Set*>,double,double,double,
            double*,double*,double*);
extern void GetOrbitalDerivatives(Spin*,int,double*,double*,double*,
            double*,double*,double*,int);
extern void GetDensityLikeDerivatives(Matrix*,vector<Set*>,double*,
            double*,double*,double*,double*,double*,int);
extern void GetSpeciesDensityDerivatives(Spin*,int,double*,double*,double*,
            double*,double*,double*,int);
extern void GetDensityDerivatives(System*,int,double*,double*,double*,double*,
            double*,double*,int);
extern void GetChargeDensityDerivatives(System*,int,double*,double*,double*,
            double*,double*,double*,int);

extern void GetShapeDerivatives(System*,int,double*,double*,double*,double*,
            double*,double*,int);
extern void GetShannonDerivatives(System*,int,double*,double*,double*,double*,
            double*,double*,int);
extern void GetElectronicSpinDensityDerivatives(System*,int,double*,double*,
            double*,double*,double*,double*,int);

#endif  // GN_EVALUATOR_H
