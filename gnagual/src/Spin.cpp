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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>

#include <Spin.h>
#include <Molecule.h>
#include <Atom.h>
#include <Element.h>
#include <Matrix.h>
#include <System.h>
#include <Units.h>
#include <Math.h>

using namespace std;

Spin::Spin(System* is,int iid, char* iname, double imass, double icharge,
int inpart)
{
  id = iid;
  npart = inpart;
  HOMO = npart-1;
  mass = imass;
  charge = icharge;
  sys = is;
  strcpy(name,iname);


  energies = 0;
  oo = 0;

  basis.clear();
};

// Build AO density matrix (P)
void Spin::BuildDensityMatrix()
{
  P[0]->SetZero();
  BuildPartialDensityMatrix(P[0],0,HOMO);
}

void Spin::BuildPartialDensityMatrix(Matrix *Q,int ll, int ul)
{
  int i,mu,nu;
  double moc[nco];

  for (i=ll;i<=ul;i++)
  {
    C->GetRowValues( i , moc );
    for (mu=0;mu<nco;mu++)
      for (nu=0;nu<nco;nu++)
        Q->ShiftValue(mu,nu,oo[i]*moc[mu]*moc[nu]);
  }
}

int Spin::Id()
{
  return id;
}
