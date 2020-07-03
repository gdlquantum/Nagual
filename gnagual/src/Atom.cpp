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
#include <iostream>

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include <Atom.h>
#include <Element.h>
#include <Units.h>
#include <Vector.h>

using namespace std;

Atom::Atom(char* il,double ix,double iy, double iz)
{
  Setup(SymbolToAtomicNumber( il ),ix,iy,iz);
};

Atom::Atom(int an,double ix,double iy, double iz)
{
  Setup(an,ix,iy,iz);
};

void Atom::Setup(int an,double ix,double iy, double iz)
{
  bonded_to.clear();

  // Position
  x = ix;
  y = iy;
  z = iz;

  // Null force
  force[0] = 0.0;
  force[1] = 0.0;
  force[2] = 0.0;

  // notice that charge is reset
  SetAtomicNumber(an);
};

void Atom::SetAtomicNumber(int an)
{
  atomic_number = an;

  // notice that charge is reset
  SetCharge( double(atomic_number) );
}


void Atom::SetSymbol(char* is)
{
  atomic_number = SymbolToAtomicNumber( is );

  // notice that charge is reset
  SetCharge( double(atomic_number) );
}

void Atom::SetCharge(double c)
{
  znuc = c;
}

double Atom::GetCovalentRadius(void)
{
  double r;

  r = ELEMENT_COV_R[atomic_number];

  return r;
};

double Atom::GetVDWRadius(void)
{
  double r;

  r = ELEMENT_VDW_R[atomic_number];

  return r;
};

static double color[3];
double* Atom::GetColor(void)
{
  //double color[3];

  color[0] = ELEMENT_COLOR[atomic_number][0];
  color[1] = ELEMENT_COLOR[atomic_number][1];
  color[2] = ELEMENT_COLOR[atomic_number][2];

  return color;
};

int Atom::SymbolToAtomicNumber( char *sym )
{
  int an = 0;
  char tst[3];

  if ( strlen(sym) == 0 )
  {
    return 0;
  }
  else
  {
    tst[0] = toupper(sym[0]);
    if ( strlen(sym) == 1 )
    {
      tst[1] = '\0';
    }
    else
    {
      if ( sym[1] >= '0' && sym[1] <= '9' )
      {
        tst[1] = '\0';
      }
      else
      {
        tst[1] = tolower(sym[1]);
        tst[2] = '\0';
      };
    };
  };

  for ( int i = 0 ; i <= GN_MAX_ATOMIC_NUMBER ; i++ )
  {
    if ( strcmp( tst , ELEMENT_SYMBOL[i] ) == 0 )
    {
      an = i;
      break;
    };
  };
  return an;
}

double Atom::Distance(Atom* other)
{
  double d;
  Vector v;

  v = Vector(this->x,this->y,this->z);
  v -= Vector(other->x,other->y,other->z);
  d = v.Norm();
  return d;
};

