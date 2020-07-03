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
#ifndef GN_SHELL_H
#define GN_SHELL_H

#include <vector>

#include <GNagual.h>
#include <Vector.h>
#include <Matrix.h>

#define PRIMARY_BASIS_SHELL      0     // Primary basis (electronic, protonic, ... any)
#define AUXILIARY_BASIS_SHELL    1     // Auxiliary basis set (electronic, ... )
#define CENTRAL_FORCE_SHELL      2     // Effective/Model/Solvent/... potentials
#define SEMILOCAL_ECP_SHELL      3     // Effective core potentials

using namespace std;

class Shell 
{
  public:

    Shell(void);

    int n;
    int l;
    int nco;
    int ll;   // Lower limit of its CAOs in basis set
    int ul;   // Upper limit of its CAOs in basis set
    int type;
    vector<double> z;
    vector<double> d;
    vector<double> ncsto;

    void Print(char*);
    void Normalize(void);
    void EvaluateRadius(void);
    //void NormalizeInteraction(Integrator*);
    double Radius(void);
    double RadialPotential(double);
    void EvaluateGTO(double*,double*,double*,double*,int,int);
    void EvaluateGTODerivative(double*,double*,double*,double*,double*,double*,int, int);
    void EvaluateGTOHessian(double*,double*,double*,double*,double*,double*,
         double*,double*,double*,int, int);

  protected:

    bool is_normalized;
    double radius;
};

#endif // GN_SHELL_H
