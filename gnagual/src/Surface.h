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
#ifndef GN_SURFACE_H
#define GN_SURFACE_H

#include <vector>

using namespace std;

class System;
class Grid;
class Set;
class Matrix;

class Surface 
{
  public:
    Surface(System*,Grid*);

    System* sys;
    Grid* grid;

    char name[256];
    int type;
    int style;
    int shiny;
    int set_number;
    int orbital_number;
    int perturbation_number;
    int point_size;
    int line_width;
    int nscalar;
    double iso;
    double color[4];
    unsigned long list;
    double* scalar;

    void Setup(void);
    void Build(void);
    void GetCubeValues(int,int,int,double*);
    void GetDerivatives(double*,double*,double*,double*,double*,double*,int);
    void GetValues(double*,double*,double*,double*,int);
    void GetExtrema(double*,double*);

  protected:

    void BuildGTO(void);
    void BuildOrbital(void);
    void BuildDensity(void);
    void BuildSpeciesDensity(void);
    void BuildChargeDensity(void);
    void BuildDensityLike(Matrix*,vector<Set*>);
    void BuildElectronicSpinDensity(void);
    void BuildShannon(void);
    void BuildShape(void);
    void BuildElectrostaticPotential(void);
};

#endif  // GN_SURFACE_H

