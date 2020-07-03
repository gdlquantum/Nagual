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
#ifndef GN_SYSTEM_H
#define GN_SYSTEM_H

#include <Parameter.h>
#include <Matrix.h>

class Molecule;
class Spin;
class QChem;

class System 
{
  public:

    System(void);

    QChem *qchem;
    Molecule* mol;
    int nspin;
    double *charges;
    double exc;
    double energy;
    double temperature;
    double chemical_potential;
    Spin* spin[GN_MAX_NSPIN];

    // Constrained SCF variables

    void AddSpin(char*,double,double,int);

    void Read(char*,char*);
    void ReadNGL(char*);
    void Setup(QChem*);
    void Print(char*);

    double Number(int,int);
    void ChangeNumber(double,int,int);
    void EvaluatePromolecularDensity(Matrix*);

    double GetAtomCharge(int);
};

#endif // GN_SYSTEM_H
