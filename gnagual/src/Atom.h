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
#ifndef GN_ATOM_H
#define GN_ATOM_H

#include <vector>

using namespace std;

class Atom 
{
  public:

    Atom(char*,double,double,double);
    Atom(int,double,double,double);

    bool basis_center;
    int atomic_number;
    int ref_bond;
    int ref_angle;
    int ref_dihedral;
    vector<int> bonded_to;  
    double x,y,z;
    double znuc;
    double charge;
    double bond;
    double angle;
    double dihedral;
    double force[3];

    void Setup(int,double,double,double);
    void SetSymbol(char*);
    void SetCharge(double);
    void SetAtomicNumber(int);
    double* GetColor(void);
    double GetCovalentRadius(void);
    double GetVDWRadius(void);
    double Distance(Atom*);

  protected:

    int SymbolToAtomicNumber( char* );
};

#endif // GN_ATOM_H
