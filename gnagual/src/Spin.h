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
#ifndef GN_SPIN_H
#define GN_SPIN_H

#include <vector>

#include <Parameter.h>
#include <Matrix.h>

using namespace std;

class System;
class Molecule;
class Set;
class File;

class Spin
{
  public:

    Spin(System*,int,char*,double,double,int);

    char name[GN_MAX_STR_SIZE];
    int id;
    int npart;
    int HOMO;
    int nco;
    double charge;
    double mass;
    double *energies;  // Orbital energies
    double *oo;        // Orbital occupation
    System* sys;           
    vector<Set*> basis;           

    Matrix *C;
    vector<Matrix*> P; // Linear response matrix

    void SetupSCF(char*,char*,int,bool);
    void Print(char*);
    void OneParticleMatrix(Matrix*,char*);
    void TwoParticleMatrix(Matrix*,Spin*,bool);
    void TwoParticleMatrixDF(Matrix*,Spin*,double*,double);
    void NewMOs(bool);
    void OrthogonalizationMatrix(void);
    void EvaluateERI4Batch(Spin*,int,int,Matrix);
    void EvaluateERI4BatchRI(Spin*,int,int,Matrix);
    void EvaluateCharges(double*,char*,bool);
    void SCFDensityMatrix(double);
    void BuildDensityMatrix(void);
    void BuildPulayMatrix(Matrix*);
    void BuildPartialDensityMatrix(Matrix*,int,int);
    void BuildFrontierFukuiMatrix(Matrix*,int,int);
    void BuildCoreMatrix(void);
    void GetDensityMatrix(Matrix*);
    void FittingVector(Matrix*,double*);
    void BuildFittingMatrix(Matrix*);
    void FittingMatrix(void);
    void FittingCoefficients(void);
    int Id(void);
    double GetOrbitalEnergy(int);
    double Entropy(void);

  protected:

    bool direct;

    void GetOccupations(void);
};

#endif // GN_SPIN_H
