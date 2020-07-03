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
#ifndef GN_MOLECULE_H
#define GN_MOLECULE_H

#include <vector>

#include <Vector.h>

using namespace std;

class Atom;

class Molecule
{
  public:
    Molecule(void);

    vector<Atom*> atom;

    void ReadXYZ(const char*);
    void ReadZMatrix(const char*);
    void ReaddeMon(bool*);

    void WriteXYZ(const char*,bool);
    void WriteZMatrix(const char*);
    void WritedeMon(void);

    void AddAtom(int,int,int,int,double,double,double);
    void DeleteLastAtom(void);
    int Natom(void);
    Atom* GetAtom(int);
    void C2Z(int); 
    void Z2C(int); 
    void Clear(void);
    void Center(void);
    void Nullify(char*);
    int PGBDrv(char*,int*); 
    int SGBDrv(char*,int*); 
    int NumberOfParticles(char*);
    double NuclearRepulsionEnergy(void);

    void BuildFromSMILES(const char*);

  protected:

    int AppIGroup(int*);
    int AppOGroup(int*);
    int AppTGroup(int*);
    int AppSGroup(int,int*);
    int AppDGroup(char*,int,int*);
    int AppCGroup(char*,int,int*);
    void AppSigma(Vector,int);
    void AppCAxis(int,int,Vector,int); 
    void AppSAxis(int,int,Vector,int);
    void AppTrans(Vector,double,int);
    void ClonAtom(int,Vector);
    double Dihedral(int,int,int,int);
};

#endif // GN_MOLECULAR_H
