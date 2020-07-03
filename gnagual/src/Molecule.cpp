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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <QtGui>

#include <Parameter.h>
#include <Molecule.h>
#include <GNagual.h>
#include <Atom.h>
#include <Units.h>
#include <Math.h>
#include <Element.h>

using std::ifstream;

Molecule::Molecule(void)
{
  atom.clear();
};

void Molecule::Clear(void)
{
  atom.clear();
}

void Molecule::ReadXYZ(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file "<<filename<<endl;
    return;
  };

  int iatom,n;
  char str[GN_MAX_STR_SIZE];
  double x,y,z;
  Atom* c;

  f >> n;
  if (n<1) cout << "WARNING in ReadXYZ: Number of atoms non positive" <<endl;
  f.getline(str,GN_MAX_STR_SIZE);
  f.getline(str,GN_MAX_STR_SIZE);
  bool inbohr;
  if (strcasestr(str,"BOHR")!=NULL) inbohr = true;
  else inbohr = false;
  cout << str << endl;
  for ( iatom = 0 ; iatom < n ; iatom++ )
  {
    f >> str;
    f >> x;
    f >> y;
    f >> z;
    if (!inbohr)
    {
      x = AngstromToBohr(x);
      y = AngstromToBohr(y);
      z = AngstromToBohr(z);
    }
    c = new Atom(str,x,y,z);
    atom.push_back( c );
  }
  Center();
 
  // Build provisional Z-Matrix
  C2Z(0);
}

int Molecule::Natom(void)
{
  return atom.size();
}

Atom* Molecule::GetAtom(int id)
{
  if (id<(signed)atom.size())
    return atom[id];
  else
  {
    exit(EXIT_FAILURE);
    return atom[0];
  }
}

void Molecule::AddAtom(int an, int r1, int r2,int r3,
                   double bond,double angle,double dihedral)
{
  Atom* a;

  a = new Atom(an,0.0,0.0,0.0);
  if (atom.size()>0)
  {
    a->ref_bond = r1;
    a->bond = bond;
    if (atom.size()>1)
    {
      a->ref_angle = r2;
      a->angle = angle;
      if (atom.size()>2)
      {
        a->ref_dihedral = r3;
        a->dihedral = dihedral;
      }
      else
      {
        a->ref_dihedral = 0;
        a->dihedral = 0.0;
      }
    }
    else
    {
      a->ref_angle = 0;
      a->angle = 0.0;
    }
    a->bonded_to.push_back( r1 );
    atom[r1]->bonded_to.push_back( atom.size() );
  } 
  else 
  {
    a->ref_bond = 0;
    a->bond = 0.0;
    a->bonded_to.clear();
  }

  atom.push_back( a );


  Z2C( atom.size()- 1 ); 
}

void Molecule::DeleteLastAtom(void)
{
  int n = atom.size();
  if (n>0) atom.resize( n-1 );
  Center();
};

void Molecule::Z2C( int llatom ) 
{
  int iatom,bondref,angleref,dihedralref;

  for ( iatom = llatom ; iatom < (signed)atom.size() ; iatom++ )
  {
    bondref = atom[iatom]->ref_bond;

    if (iatom>0)
    {
      atom[iatom]->x = atom[bondref]->x;
      atom[iatom]->y = atom[bondref]->y;
      atom[iatom]->z = atom[bondref]->z;
    }
    else
    {
      atom[iatom]->x = 0.0;
      atom[iatom]->y = 0.0;
      atom[iatom]->z = 0.0;
    }

    if ( iatom > 0 )
    {
      Vector shift(0.0,0.0,1.0);
      if ( iatom > 1 )
      {
        angleref = atom[iatom]->ref_angle;
        shift.x = atom[angleref]->x - atom[bondref]->x;
        shift.y = atom[angleref]->y - atom[bondref]->y;
        shift.z = atom[angleref]->z - atom[bondref]->z;
        shift ^= (1.0);    
      };
      shift *= atom[iatom]->bond;
      atom[iatom]->x += shift.x;
      atom[iatom]->y += shift.y;
      atom[iatom]->z += shift.z;

      if ( iatom > 1 )
      {
        Vector axis;
        shift.x = atom[angleref]->x - atom[bondref]->x;
        shift.y = atom[angleref]->y - atom[bondref]->y;
        shift.z = atom[angleref]->z - atom[bondref]->z;
        if ( iatom == 2 )
        {
          axis = Vector(0.0,1.0,0.0) < shift;
        }
        else
        {
          dihedralref = atom[iatom]->ref_dihedral;
          axis.x = atom[dihedralref]->x - atom[bondref]->x;
          axis.y = atom[dihedralref]->y - atom[bondref]->y;
          axis.z = atom[dihedralref]->z - atom[bondref]->z;
          axis <= shift;
        };
        shift.x = atom[iatom]->x - atom[bondref]->x;
        shift.y = atom[iatom]->y - atom[bondref]->y;
        shift.z = atom[iatom]->z - atom[bondref]->z;
        if (shift.Colineal(axis)) axis = shift.Orthogonal();
        shift.Rotate(axis,atom[iatom]->angle);
        atom[iatom]->x = atom[bondref]->x + shift.x;
        atom[iatom]->y = atom[bondref]->y + shift.y;
        atom[iatom]->z = atom[bondref]->z + shift.z;

        if ( iatom > 2 )
        {
          axis.x = atom[angleref]->x - atom[bondref]->x;
          axis.y = atom[angleref]->y - atom[bondref]->y;
          axis.z = atom[angleref]->z - atom[bondref]->z;
          shift.x = atom[iatom]->x - atom[bondref]->x;
          shift.y = atom[iatom]->y - atom[bondref]->y;
          shift.z = atom[iatom]->z - atom[bondref]->z;
          shift.Rotate(axis,atom[iatom]->dihedral);
          atom[iatom]->x = atom[bondref]->x + shift.x;
          atom[iatom]->y = atom[bondref]->y + shift.y;
          atom[iatom]->z = atom[bondref]->z + shift.z;
        };
      };
    };
  };
  Center();
};


void Molecule::Center(void)
{
    
  Vector ref(0.0,0.0,0.0);
  for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    ref.x += atom[iatom]->x;
    ref.y += atom[iatom]->y;
    ref.z += atom[iatom]->z;
  }
  ref /= (double)(atom.size());
  for (int iatom = 0 ; iatom < (signed) atom.size() ; iatom++ )
  {
    atom[iatom]->x -= ref.x;
    atom[iatom]->y -= ref.y;
    atom[iatom]->z -= ref.z;
  }
  
};

// Point Group operations applier driver
int Molecule::PGBDrv(char* pg,int* ref) 
{
  int n;

  if ( strcmp(pg,"Cs") == 0 )
  {
    strcpy(pg,"C1h");
  }
  else if ( strcmp(pg,"Ci") == 0 )
  {
    strcpy(pg,"S2");
  }
  else if ( strcmp(pg,"Vh") == 0 )
  {
    strcpy(pg,"D2h");
  } 
  else if ( strcmp(pg,"Vd") == 0 )
  {
    strcpy(pg,"D2d");
  }
  else if ( strcmp(pg,"S6v") == 0 )
  {
    strcpy(pg,"D3d");
  }
  else if ( strcmp(pg,"S8v") == 0 )
  {
    strcpy(pg,"D4d");
  };

  n = 0;
  if ((pg[0] != 'T')&&(pg[0] != 'O')&&(pg[0] != 'I')&&
      (pg[1] != '*')&&(pg[1] != 's'))
  {
    char s[10];
    sprintf(s,"%c",pg[1]);
    n = atoi(s);
  };

  Vector v(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  for ( int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    atom[iatom]->x -= v.x;
    atom[iatom]->y -= v.y;
    atom[iatom]->z -= v.z;
  }

  if ( pg[0] == 'C' )
  {
    AppCGroup(pg,n,&ref[1]);
  }
  else if ( pg[0] == 'D' )
  {
    AppDGroup(pg,n,&ref[1]);
  }
  else if ( pg[0] == 'S' )
  {
    AppSGroup(n,&ref[1]);
  }
  else if ( pg[0] == 'T' )
  {
    AppTGroup(&ref[1]);
  }
  else if ( pg[0] == 'O' )
  {
    AppOGroup(&ref[1]);
  }
  else if ( pg[0] == 'I' )
  {
    AppIGroup(&ref[1]);
  };

  Center();

  return 0;
};


// Purpose: Space Group operations applier driver
int Molecule::SGBDrv(char* sg,int* ref) 
{
  int nr;

  nr = atom.size();

  if ( sg[1] == '1' )
  {
    Vector axis(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
    axis.x -= atom[ref[0]]->x;
    axis.y -= atom[ref[0]]->y;
    axis.z -= atom[ref[0]]->z;
    double sd = QInputDialog::getDouble(0,
                "Translation Unit Length","R = ",1.0,0.0,10000.0,3);
    double nt = QInputDialog::getInt(0,
                "Number on Monomers","NM = ",2,2,10000);
    nt = nt - 1;
    double s = 0.0;
    for ( int it = 0 ; it < nt ; it++ )
    {
      s = s + sd;
      AppTrans(axis,s,nr);
    }
  }
  else 
  {
    cout << "sgbdrv: Unsupported spatial group"<<endl;
    exit(EXIT_FAILURE);
  };

  Center();
  return 0;
};

// Apply I group operations 
int Molecule::AppIGroup(int *ref)
{
  int i,iaxis,nr;

  nr = atom.size();

  Vector v,axisa,axisb,axis[12];
  Vector axis2[20],axis3[15],plane[15];

  axisa = Vector(atom[ref[0]]->x,atom[ref[0]]->y, atom[ref[0]]->z);
  axisa /= axisa.Norm();
  axisb = Vector(atom[ref[1]]->x,atom[ref[1]]->y, atom[ref[1]]->z);
  v = Vector( axisa );
  v *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= v;
  axisb /= axisb.Norm();
  v = axisa > axisb;

  axis[0] = axisa;
  axis[1] = axisa;
  axis[1].Rotate( v , atan(2.0)*180.0/M_PI );
  for ( iaxis = 2 ; iaxis < 6 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 72.0 );
  };
  axis[6] = axis[4];
  axis[6] *= -1.0; 
  axis[7] =  axis[5];
  axis[7] *= -1.0; 
  axis[8] = axis[1];
  axis[8] *= -1.0; 
  axis[9] = axis[2];
  axis[9] *= -1.0; 
  axis[10] = axis[3];
  axis[10] *= -1.0; 
  axis[11] = axis[0];
  axis[11] *= -1.0; 

  for ( iaxis = 0 ; iaxis < 12 ; iaxis++ )
  {
    for ( i = 0 ; i < 4 ; i++ )
    {
      AppSAxis( 10 , i + 1 , axis[iaxis] , nr );
    };
  };

  axis2[0] = axis[0] + axis[1] + axis[2];
  axis2[0] /= axis2[0].Norm();
  axis2[5] = axis[1] + axis[2] + axis[6];
  axis2[5] /= axis2[5].Norm();
  for ( iaxis = 1 ; iaxis < 5 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-1];
    axis2[iaxis].Rotate( axisa , 72.0 );
    axis2[iaxis+5] = axis2[iaxis+4];
    axis2[iaxis+5].Rotate( axisa , 72.0 );
  };
  for ( iaxis = 10 ; iaxis < 20 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-10];
    axis2[iaxis] *= -1.0; 
  };

  for ( iaxis = 0 ; iaxis < 20 ; iaxis++ )
  {
    for ( i = 0 ; i < 2 ; i++ ) 
    {
      AppSAxis( 6 , i + 1 , axis2[iaxis] , nr );
    };
  };

  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    axis3[iaxis] = axis[0] + axis[iaxis+1];
    axis3[iaxis+5] = axis[iaxis+1] + axis[iaxis+2];
  };
  axis3[4] = axis[0] + axis[5];
  axis3[9] = axis[5] + axis[1];
  axis3[10] = axis[1] + axis[6];
  for ( iaxis = 11 ; iaxis < 15 ; iaxis++ )
  {
    axis3[iaxis] = axis3[iaxis-1];
    axis3[iaxis].Rotate( axisa , 72.0 ); 
  };

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    AppCAxis( 2 , 1 , axis3[iaxis] , nr );
  };

  for (iaxis = 0 ; iaxis < 5 ; iaxis++ )
  {
    plane[iaxis] = axisa > axis[iaxis+1];
    plane[iaxis+10] = axis[iaxis+1] > axis[iaxis+6];
  };
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    plane[iaxis+5] = axis[iaxis+1] > axis[iaxis+2];
  };
  plane[9] = axis[5] > axis[1];

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    AppSigma( plane[iaxis] , nr );
  };

  AppSAxis( 2 , 1 , axisa , nr );

  return 0;
};


// Apply O group operations 
int Molecule::AppOGroup(int *ref)
{
  int i,j,iaxis;    
  int nr;            
  Vector axisa;    
  Vector c2axis;  
  Vector axis[5];

  nr = atom.size();

  axis[0] = Vector(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axis[0] /= axis[0].Norm();
  axisa = axis[0];
  axis[1] = Vector(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  axis[2] = Vector( axis[0] );
  axis[2] *= axis[1].Dot( axis[0] )/axis[0].Dot( axis[0] );
  axis[1] -= axis[2];
  axis[1] /= axis[1].Norm();
  axis[2] = axis[0] > axis[1];
  axis[3] = axis[0];
  axis[4] = axis[1];

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    c2axis = axis[iaxis+1] + axis[iaxis+2];
    c2axis /= c2axis.Norm(); 
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSAxis( 4 , i+1, axis[iaxis] , nr );
      AppCAxis( 4 , i+1, axis[iaxis] , nr );
      AppCAxis( 2 , 1, c2axis , nr );
      c2axis.Rotate( axis[iaxis] , 90.0 );
      AppCAxis( 2 , 1, c2axis , nr );
      for ( j = 0 ;  j < 3 ; j++ )
      {
        AppSigma( axis[iaxis+j] , nr );
      };
    };
  };

  AppSAxis( 2 , 1 , axisa , nr );

  axis[0] += axis[1];
  axis[0] += axis[2];
  axis[0] /= axis[0].Norm();
  for ( iaxis = 1 ; iaxis < 4 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 90.0 );
  };
 
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    for ( i = 0 ; i < 5 ; i++ )
    {
      AppSAxis( 6 , i+1, axis[iaxis] , nr );
    };
  };

  return 0;
};

// Apply T group operations 
int Molecule::AppTGroup(int *ref)
{
  int       i,iaxis,nr;
  double    tangle; 
  Vector axisb;  
  Vector v;
  Vector axis[4]; 

  nr = atom.size();

  axis[0] = Vector(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axis[0] /= axis[0].Norm();
  axisb = Vector(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  axis[2] = Vector( axis[0] );
  axis[2] *= axisb.Dot( axis[0] )/axis[0].Dot( axis[0] );
  axisb -= axis[2];
  axisb /= axisb.Norm();

  tangle = acos(-1.0/3.0);
  axis[2] = axis[0] > axisb;
  axis[1] = axis[0];
  axis[1].Rotate( axis[2] , 180.0*tangle/M_PI );
  axis[2] = axis[1];
  axis[2].Rotate( axis[0] , 120.0 );
  axis[3] = axis[2];
  axis[3].Rotate( axis[0] , 120.0 );

  for ( iaxis = 0 ; iaxis  < 4 ; iaxis++ )
  {
    AppCAxis(3,1, axis[iaxis], nr );
    AppCAxis(3,2, axis[iaxis], nr );
    if ( iaxis == 0 )
    {
      v = axis[0] > axisb;
    }
    else 
    {
      v = axis[iaxis-1] > axis[iaxis];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSigma( v , nr );
      v.Rotate( axis[iaxis] , 120.0 );
    };
  };

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    if ( iaxis == 2 )
    {
      v = axis[iaxis] + axis[1];
    }
    else
    {
      v = axis[iaxis] + axis[iaxis+1];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      AppSAxis(4, i + 1, v , nr );
    };
  };

  return 0;
};

// Apply S group operations 
int Molecule::AppSGroup(int n,int *ref)
{
  int  i;    
  int nr;    
  char s[10];
  
  if ( n%2 == 1 )
  {
    sprintf(s,"C%dh",n);
    return AppCGroup(s,n,ref);
  }
  else
  {
    nr = atom.size();

    Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
    axisa /= axisa.Norm();

    if ( n >= 4 )
    {
      for ( i = 1 ; i <= n/2-1 ; i++ )
      {
        AppCAxis( n , i , axisa , nr );
      };
    };

    for ( i = 1 ; i < n ; i++ )
    {
      AppSAxis( n , i , axisa , nr );
    };

  };

  return 0;
};

// Apply D group operations 
int Molecule::AppDGroup(char* pg,int n,int *ref)
{
  int i,top,step,nr; 

  nr = atom.size();

  Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);
  axisa /= axisa.Norm();

  if ( pg == (char*)"D*h" )
  {
    AppSAxis(2,1,axisa,nr);
    return 0;
  };
  
  for ( i = 1 ; i < n ; i++ )
  {
    AppCAxis( n , i , axisa , nr );
  };

  Vector axisb(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
  Vector v( axisa );
  v *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= v;
  axisb /= axisb.Norm();
  v = axisb;
  for ( i = 1 ; i <= n ; i++ )
  {
    AppCAxis( 2 , 1 , v , nr );
    v.Rotate( axisa , 180.0/((double)n) );
  };

  if ( pg[2] == 'h' )
  {
    if ( n%2 == 0 )
    {
      top = n - 1;
      step = 1;
    }
    else
    {
      top = 2*n-1;
      step = 2;
    };
    for ( i = 1 ; i <= top ; i += step )
    {
      AppSAxis( n , i , axisa , nr );
    };
  }
  else if ( pg[2] == 'd' )
  {
    for ( i = 1 ; i <= 2*n-1 ; i += 2 )
    {
      AppSAxis( n , i , axisa , nr );
    };
  };

  v = axisa > axisb;
  if ( pg[2] == 'd' )
  {
    v.Rotate( axisa , 180.0/((double)n) );
  };
  for ( i = 1 ; i <= n ; i++ )
  {
    AppSigma( v , nr );
    v.Rotate( axisa , 360.0/((double)n) );
  };

  if ( pg[2] == 'h' )
    AppSigma( axisa , nr );

  return 0;
};

// Apply C group operations 
int Molecule::AppCGroup(char* pg,int n,int *ref)
{
  int i,nr;  

  if ( n < 2 ) return 0;

  Vector axisa(atom[ref[0]]->x,atom[ref[0]]->y,atom[ref[0]]->z);

  nr = atom.size();
  for ( i = 1 ; i < n ; i++ )
  {
    AppCAxis(n,i,axisa,nr);
  };
  
  if ( strlen(pg) == 2 ) return 0;
  
  if ( pg[2] == 'h' )
  {
    for ( i = 1 ; i < n ; i++ )
    {
      AppSAxis( n , i , axisa , nr );
    };
    AppSigma( axisa , nr );
  }
  else 
  {
    Vector axisb(atom[ref[1]]->x,atom[ref[1]]->y,atom[ref[1]]->z);
    Vector normal( axisa );
    normal *= axisb.Dot( axisa )/axisa.Dot( axisa );
    axisb -= normal;
    normal = axisa > axisb;
    for ( i = 1 ; i <= n ; i++ )
    {
      AppSigma(normal,nr);
      normal.Rotate(axisa,360.0/((double)n));
    };
  };
  return 0;
};


// Apply sigma plane reflexion 
void Molecule::AppSigma(Vector normal,int nr)
{
  int iatom;
  double sp;
  Vector va,vb;

  normal ^= 1.0;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    va = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    sp = va.Dot(normal);
    vb = normal;
    vb *= 2.0*sp;
    va  -= vb;
    ClonAtom(iatom,va);
  };
};

// Apply Cn rotation operation 
void Molecule::AppCAxis(int n,int i,Vector axis,int nr) 
{
  int iatom; 
  double phi; 
  Vector v;

  phi = 360.0*((double)i)/((double)n);
  axis ^= 1.0;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    v.Rotate(axis,phi); 

    ClonAtom(iatom,v);
  };
};

// Apply Sn rotation operation 
void Molecule::AppSAxis(int n,int i,Vector axis,int nr)
{
  int iatom;
  double sp,phi;
  Vector v,va,vb;

  if ( i%2 == 0 )
  {
    AppCAxis(n,i,axis,nr);
  }
  else
  {
    phi = 360.0*((double)i)/((double)n);
    axis ^= 1.0;

    for ( iatom = 0 ; iatom < nr ; iatom++ )
    {
      if ( atom[iatom]->atomic_number == 0 ) continue;

      v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
      v.Rotate(axis,phi); 
      sp = Vector( v ).Dot( axis );
      va = Vector( v );
      vb = axis;
      vb *= 2.0*sp;
      va  -= vb;
      ClonAtom(iatom,va);
    };
  };
}; 

// Apply translation
void Molecule::AppTrans(Vector axis,double shift, int nr)
{
  int iatom;
  Vector normal,v;

  normal = axis;
  normal ^= 1.0;
  normal *= shift;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;

    v = Vector(atom[iatom]->x,atom[iatom]->y,atom[iatom]->z);
    v += normal;
    ClonAtom(iatom,v);
  };
};

// Make a copy of a non-redundant atom if required
void Molecule::ClonAtom(int from,Vector newpos)
{
  Vector v;  
  
  for ( int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    if ( atom[iatom]->atomic_number == 0 ) continue;
    v = newpos;
    v.x -= atom[iatom]->x;
    v.y -= atom[iatom]->y;
    v.z -= atom[iatom]->z;
    if ( v.Norm() < 0.5 )
    {
      return;
    };
  };

  Atom *a = new Atom(atom[from]->atomic_number,newpos.x,newpos.y,newpos.z);
  atom.push_back( a );
};

void Molecule::WriteXYZ(const char *filename, bool head)
{
  ofstream f;

  if (head)
  {
    f.open(filename);
    f << atom.size() <<endl;
    f << endl;
  }
  else f.open(filename,ios::app);

  for (int i=0;i<(signed)atom.size();i++)
  {
    f << ELEMENT_SYMBOL[atom[i]->atomic_number] 
         << fixed << setw(20) << setprecision(10) << "  " 
      << BohrToAngstrom(atom[i]->x) << "   "
      << BohrToAngstrom(atom[i]->y) << "   "
      << BohrToAngstrom(atom[i]->z) << endl;
  }
  f << endl;
  f.close();
};

int Molecule::NumberOfParticles(char* name)
{
  if ( strncasecmp(name,"ELECT",5) == 0 ) 
  {
    int nelec;
  
    nelec = 0;
    for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
    {
      nelec += atom[iatom]->atomic_number;
    }
    return nelec;
  }
  else if ( strncasecmp(name,"PROTO",5) == 0 ) 
  {
    int nh;
  
    nh = 0;
    for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
    {
      if (atom[iatom]->atomic_number==1) nh++;
    }
    return nh;
  } 
  else
  {
    return 0;
  }
};

// Comput nuclear repulsion energy in atomic units
double Molecule::NuclearRepulsionEnergy()
{
  int iatom,jatom;
  double d,nr,x,y,z;

  nr = 0.0;
  for (iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
  {
    for (jatom = iatom+1 ; jatom < (signed)atom.size() ; jatom++ )
    {
      if (atom[iatom]->znuc>0&&atom[jatom]->znuc>0)
      {
        x = atom[iatom]->x - atom[jatom]->x;
        y = atom[iatom]->y - atom[jatom]->y;
        z = atom[iatom]->z - atom[jatom]->z;
        d = sqrt(x*x+y*y+z*z);
        nr += atom[iatom]->znuc*atom[jatom]->znuc/d;
      }
    }
  }
  return nr;
}

void Molecule::Nullify(char* name)
{
  if ( strncasecmp(name,"PROTO",5) == 0 ) 
  {
    for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
    {
      if (atom[iatom]->atomic_number==1) 
        atom[iatom]->znuc = 0.0;
    }
  }
  else if ( strncasecmp(name,"SOLVE",5) == 0 ) 
  {
    for (int iatom = 0 ; iatom < (signed)atom.size() ; iatom++ )
    {
      if (atom[iatom]->atomic_number==104) 
        atom[iatom]->znuc = 0.0;
    }
  }
};

void Molecule::C2Z( int llatom ) 
{
  int iatom;
  double val;
  Vector va,vb;

  for ( iatom = llatom ; iatom < (signed)atom.size() ; iatom++ )
  {
    atom[iatom]->ref_bond = 0;
    atom[iatom]->bond = 0.0;
    atom[iatom]->ref_angle = 0;
    atom[iatom]->angle = 0.0;
    atom[iatom]->ref_dihedral = 0;
    atom[iatom]->dihedral = 0.0;

    if (iatom>0)
    {
      va = Vector(atom[iatom]->x - atom[iatom-1]->x,
                  atom[iatom]->y - atom[iatom-1]->y,
                  atom[iatom]->z - atom[iatom-1]->z);
      atom[iatom]->ref_bond = iatom - 1;
      atom[iatom]->bond = va.Norm();
      if (iatom>1)
      {
        va ^= 1.0;
        vb = Vector(atom[iatom-2]->x - atom[iatom-1]->x,
                    atom[iatom-2]->y - atom[iatom-1]->y,
                    atom[iatom-2]->z - atom[iatom-1]->z);
        vb ^= 1.0;
        atom[iatom]->ref_angle = iatom - 2;
        atom[iatom]->angle = RadianToDegree(acos(va.Dot(vb)));
        if (iatom>2) 
        {
          atom[iatom]->ref_dihedral = iatom - 3;
          val = Dihedral(iatom-3,iatom-2,iatom-1,iatom);
          atom[iatom]->dihedral = val;
        }
      }
    }
    if (atom[iatom]->angle==180.0) 
    {
      va = Vector(atom[iatom]->x - atom[iatom-1]->x,
                  atom[iatom]->y - atom[iatom-1]->y,
                  atom[iatom]->z - atom[iatom-1]->z);
      vb = va.Orthogonal();
      vb *= va.Norm();
      va = Vector(atom[iatom-1]->x,atom[iatom-1]->y,atom[iatom-1]->z);
      va += vb;
      std::vector<Atom*>::iterator it;
      it = atom.begin();
      it = atom.insert(it+iatom,new Atom(0,va.x,va.y,va.z));
      iatom--;
    }
  }
}
  
double Molecule::Dihedral(int iatom, int jatom, int katom, int latom)
{
  double val;
  double ra[3],rb[3],rc[3],rd[3];
  Vector auxa,auxb,auxc,va,vb,vc,vd,rva,rvc,rvd;

  ra[0] = atom[iatom]->x;
  ra[1] = atom[iatom]->y;
  ra[2] = atom[iatom]->z;
  va = Vector(ra);

  rb[0] = atom[jatom]->x;
  rb[1] = atom[jatom]->y;
  rb[2] = atom[jatom]->z;
  vb = Vector(rb);

  rc[0] = atom[katom]->x;
  rc[1] = atom[katom]->y;
  rc[2] = atom[katom]->z;
  vc = Vector(rc);

  rd[0] = atom[latom]->x;
  rd[1] = atom[latom]->y;
  rd[2] = atom[latom]->z;
  vd = Vector(rd);

  rva = va - vb;
  rvc = vc - vb;
  rvd = vd - vc;
  rva ^= 1.0;
  rvc ^= 1.0;
  rvd ^= 1.0;

  if (rva.Colineal(rvc)) 
  {
    if (rvd.Dot(rvc)<0.0) return 0.0; 
    else return 180.0; 
  }

  auxa = rvc;
  auxa *= -rva.Dot( rvc );
  auxa += rva;
  auxa ^= 1.0;

  auxb = rvc;
  auxb *= -rvd.Dot( rvc );
  auxb += rvd;
  auxb ^= 1.0;

  auxc = auxa;
  auxc.Rotate( rvc ,-90.0 );
  val = auxa.Dot( auxb );

  // Avoid running out of domain of acos()
  val = double(int(100000.0*val+0.5))/100000.0; 

  if ( auxb.Dot( auxc ) < 0.0 )
  {
    val = 360.0 - 180.0/M_PI*acos( val );
  }
  else
  {
    val =  180.0/M_PI*acos( val );
  }
  return val;
}

void Molecule::WriteZMatrix(const char *filename)
{
  ofstream f(filename);
  f << atom.size() <<endl;
  f << endl;
  for (int i=0;i<(signed)atom.size();i++)
  {
    f << ELEMENT_SYMBOL[atom[i]->atomic_number] << "  "; 
    if (i>0) 
      f << fixed << setw(4) << atom[i]->ref_bond << "  " 
        << fixed << setw(10) << setprecision(3) 
        << BohrToAngstrom(atom[i]->bond) << "   ";
    if (i>1) 
      f << fixed << setw(4) << atom[i]->ref_angle << "  "
        << fixed << setw(10) << setprecision(2) << atom[i]->angle << "   ";
    if (i>2) 
      f << fixed << setw(4) << atom[i]->ref_dihedral << "  "
        << fixed << setw(10) << setprecision(2) << atom[i]->dihedral;
    f << endl;
  }
  f.close();
};

void Molecule::ReadZMatrix(const char *filename)
{
  ifstream f(filename);
  if ( f.fail() ) return;

  char line[GN_MAX_STR_SIZE];
  char str[GN_MAX_STR_SIZE];
  int iatom,n;
  Atom* c;

  atom.clear();
  f >> n;
  f.getline(line,GN_MAX_STR_SIZE);
  f.getline(line,GN_MAX_STR_SIZE);
  for ( iatom = 0 ; iatom < n ; iatom++ )
  {
    f >> str;
    atom.push_back( new Atom(str,0.0,0.0,0.0) );
    c = atom[iatom];
    if (iatom>0) 
    {
      f >> c->ref_bond;
      f >> c->bond;
      c->bond = AngstromToBohr(c->bond);
    }
    if (iatom>1) 
    {
      f >> c->ref_angle;
      f >> c->angle;
    }
    if (iatom>2) 
    {
      f >> c->ref_dihedral;
      f >> c->dihedral;
    }
  }
  f.close();

  // Obtain cartesian representation
  Z2C(0); 
}

void Molecule::ReaddeMon(bool *optimized)
{
  int natom;
  double x,y,z;
  char str[GN_MAX_STR_SIZE];
  Atom* c;

  ifstream f("deMon.out");
  if ( f.fail() ) return;
  while (true)
  {
    f.getline(str,GN_MAX_STR_SIZE);
    if ( strncasecmp(str," NUMBER OF ATOMS: ",18) == 0 ) 
    {
      natom = atoi(&str[18]);
    }
    if ( strncasecmp(str," *** CONVERGED ",15) == 0 ) 
    {
      *optimized = true;
    }
    if ( ( strncasecmp(str," COORDINATES OF OPTIMIZATION STEP 1 ",36) == 0 ) 
       ||( strncasecmp(str," FINAL INPUT ORIENTATION ",25) == 0 ) ) 
    {
      f >> str; 
      f.getline(str,GN_MAX_STR_SIZE);
      atom.clear();
      for (int iatom=0;iatom<natom;iatom++)
      {
        f >> str; 
        f >> str; 
        f >> x;
        f >> y;
        f >> z;
        x = AngstromToBohr(x);
        y = AngstromToBohr(y);
        z = AngstromToBohr(z);
        c = new Atom(str,x,y,z);
        atom.push_back( c );
        f.getline(str,GN_MAX_STR_SIZE);
      }
      break;
    }
  }
  f.close();

  Center();
 
  // Build provisional Z-Matrix
  C2Z(0);
}

// Parece funcionar para CHON
// Roberto Flores-Moreno, 10-10-2017 (Utrech)
void Molecule::BuildFromSMILES(const char* smiles)
{
  char str[GN_MAX_STR_SIZE];
  char symbol[3];
  int iatom = -1;
  vector<int> dconst;
  vector<int> branch_point;
  vector<int> cycle_point;
  vector<int> cycle_point_label;
  vector<int> ring;
  Atom* c;
  int bond_ref;
  int angle_ref;
  int dihedral_ref;

  double bond;// = AngstromToBohr(1.0);
  double angle;
  double dihedral;

  bond_ref = -1;
  ring.clear();
  dconst.clear();

  int doring = 0;
  int order = 1;
  atom.clear();
  branch_point.clear();
  cycle_point.clear();
  cycle_point_label.clear();
  dihedral = 60.0;
  int i=0;
  while (smiles[i]!='\0')
  {
    if (smiles[i]=='(') // Save last branch point
    {
      branch_point.push_back( bond_ref );
      dihedral = 180.0;
      doring += 1;
    }
    else if (smiles[i]==')') // Recover to branch point
    {
      bond_ref = branch_point[branch_point.size()-1];
      dihedral += 120.0;
      branch_point.pop_back( );
      if (ring.size()>0) doring -= 1;
    }
    else if (smiles[i]=='=') 
    {
      order = 2;
    }
    else if (smiles[i]=='#') 
    {
      order = 3;
    }
    else if ((smiles[i]>='A')&&(smiles[i]<='Z'))
    {
      int j=0;
      symbol[j] = smiles[i];
      j++;
      if ((smiles[i+1]>='a')&&(smiles[i+1]<='z'))
      {
        symbol[j] = smiles[i+1];
        j++;
      }
      symbol[j] = '\0';
      c = new Atom(symbol,0.0,0.0,0.0);
      if (atom.size()>0)
      {
        bond = ELEMENT_COV_R[c->atomic_number] + 
               ELEMENT_COV_R[atom[bond_ref]->atomic_number];
      }
      angle = 109.467;
      if (atom.size()==0)
        AddAtom(c->atomic_number,-1,-1,-1,0.0,0.0,0.0);
      else if (atom.size()==1)
        AddAtom(c->atomic_number,bond_ref,-1,-1,bond,0.0,0.0);
      else if (atom.size()==2)
        AddAtom(c->atomic_number,bond_ref,1-bond_ref,-1,bond,angle,0.0);
      else 
      {
        angle_ref = atom[bond_ref]->bonded_to[0];
        if (atom[angle_ref]->bonded_to[0]==bond_ref)
        {
          dihedral_ref = atom[angle_ref]->bonded_to[1];
        }
        else
          dihedral_ref = atom[angle_ref]->bonded_to[0];
        AddAtom(c->atomic_number,bond_ref,angle_ref,dihedral_ref,
            bond,angle,dihedral);
      }
      for (size_t iorder=1;iorder<order;iorder++)
      {
        atom[atom.size()-1]->bonded_to.push_back( bond_ref );
        atom[bond_ref]->bonded_to.push_back( atom.size()-1 );
        atom[atom.size()-1]->bond -= 0.2;
      }
      order = 1;
      bond_ref = atom.size()-1;
      if ((doring==0)&&(ring.size()>0)) 
      {
        ring.push_back( bond_ref );
        if (ring.size()>3) dconst.push_back( true );
        else dconst.push_back( false );
      }
      else dconst.push_back( false );
    }
    else if ((smiles[i]>='0')&&(smiles[i]<='9'))
    {
      // Get the number
      int first = i;
      int last = first+1;
      str[0] = smiles[first];
      while ((smiles[last]>='0')&&(smiles[last]<='9')) 
      {
        str[last-first] = smiles[last];
        last++;
      }
      str[last] = '\0';
      int label = atoi(str);
     
      // Check if it was listed
      int pos = -1;
      for (size_t il=0;il<cycle_point_label.size();il++)
      {
        if (cycle_point_label[il]==label) pos = il;
      }

      // Opening ring
      if (pos<0)
      {
        cycle_point.push_back( int(atom.size()-1) );
        cycle_point_label.push_back( label );
        ring.clear();
        ring.push_back( atom.size() - 1 );
        doring = 0;
      }
      // Closing ring
      else
      {
        int cycle_ref = cycle_point[pos];
        int n = ring.size();
        cycle_point.pop_back( );
        cycle_point_label.pop_back( );
        // Build the ring by correcting previous angles
        for (size_t ir=1;ir<ring.size();ir++)
        {
          cout << "Atom in ring " << ring[ir] << endl;
          atom[ring[ir]]->angle = 180.0 - 360.0/n;
          if (ir>2) atom[ring[ir]]->dihedral = 0.0;
        }
        atom[cycle_ref]->bonded_to.push_back( atom.size()-1 );
        atom[atom.size()-1]->bonded_to.push_back( cycle_ref );
        ring.clear();
        doring = 1;
      }
    }
    i++;
  }


  // Re-conform to avoid ugly arrangement
  for (size_t jatom=3;jatom<atom.size();jatom++)
  {
    for (size_t iatom=0;iatom<jatom;iatom++)
    {
      if (atom[jatom]->Distance(atom[iatom])<1.0) 
      {
        if (!dconst[jatom])
        {
          atom[jatom]->dihedral += 120.0;
          if (atom[jatom]->dihedral>360) atom[jatom]->dihedral -= 360;
        }
      }
    }
  }

  // Add hydrogens (solo funciona para C,N,O)
  cout << "Adjusting coordination" <<endl;
  for (size_t iatom=0;iatom<atom.size();iatom++)
  {
    bond = ELEMENT_COV_R[atom[iatom]->atomic_number] + ELEMENT_COV_R[1];
    int cn = atom[iatom]->bonded_to.size();
    //cn += atom[iatom]->atomic_number-6;
    if ((atom[iatom]->atomic_number>=6)&&
        (atom[iatom]->atomic_number<=8)&&
        (cn!=4))
    {
      cout << iatom << "  " << atom[iatom]->bonded_to.size() << endl;
      cout << "CN =  " <<cn << endl;
      int other1,other2;
      bond_ref = iatom;
      symbol[0] = 'H';
      symbol[1] = '\0';
      if (cn>4) cout << "Coordination number is to high!" << endl;
      if (cn==0) 
      {
        cout << "I hope this is methane!"<<endl;
        AddAtom(1,bond_ref,1,2,bond,109.467,0.0);
        AddAtom(1,bond_ref,1,2,bond,109.467,0.0);
        if (atom[iatom]->atomic_number<8)
        {
          AddAtom(1,bond_ref,1,2,bond,109.467,120.0);
          if (atom[iatom]->atomic_number<7)
            AddAtom(1,bond_ref,1,2,bond,109.467,240.0);
        }
      } 
      else if (cn==1)
      {
        other1 = atom[iatom]->bonded_to[0];
        if (atom[other1]->bonded_to.size()>1)
        {
          if (iatom==atom[other1]->bonded_to[0])
            AddAtom(1,bond_ref,other1,atom[other1]->bonded_to[1],bond,109.467,0.0);
          else
            AddAtom(1,bond_ref,other1,atom[other1]->bonded_to[0],bond,109.467,0.0);
        }
        else AddAtom(1,bond_ref,other1,2,bond,109.467,0.0);
        other2 = atom.size() - 1;
        if (atom[iatom]->atomic_number<8)
        {
          AddAtom(1,bond_ref,other1,other2,bond,109.467,120.0);
          if (atom[iatom]->atomic_number<7)
            AddAtom(1,bond_ref,other1,other2,bond,109.467,240.0);
        }
      } 
      else if (cn==2)
      {
        if (atom[iatom]->atomic_number<8)
        {
          other1 = atom[iatom]->bonded_to[0];
          other2 = atom[iatom]->bonded_to[1];
          if (other1!=other2)
          {
            AddAtom(1,iatom,other1,other2,bond,109.467,120.0);
            if (atom[iatom]->atomic_number<7)
              AddAtom(1,iatom,other1,other2,bond,109.467,240.0);
/*
                    atom[other2]->ref_bond,
                    atom[other2]->ref_angle,
                    atom[other2]->ref_dihedral,
                    bond,109.467,atom[other2]->dihedral+120.0);
            if (atom[iatom]->atomic_number<7)
              AddAtom(1,
                    atom[other2]->ref_bond,
                    atom[other2]->ref_angle,
                    atom[other2]->ref_dihedral,
                    bond,109.467,atom[other2]->dihedral+240.0);
*/
          } 
          else 
          {
            AddAtom(1,bond_ref,other1,other2,bond,120.0,180.0);
            if (atom[iatom]->atomic_number<7)
              AddAtom(1,bond_ref,other1,other2,bond,120.0,0.0);
          }
        }
      } 
      else if (atom[iatom]->atomic_number<7)
      {
        other1 = atom[iatom]->bonded_to[0];
        other2 = atom[iatom]->bonded_to[1];
        int other3 = atom[iatom]->bonded_to[2];
        if ((other1!=other2)&&(other1!=other3)&&(other2!=other3))
        {
          AddAtom(1,bond_ref,other1,other2,bond,109.467,240.0);
        } 
        else 
        {
          if ((other3==other1)&&(other3==other2)&&(other1==other2))
            AddAtom(1,bond_ref,other1,other2,bond,180.0,180.0);
          else if (other1==other2)
            AddAtom(1,bond_ref,other1,other3,bond,120.0,180.0);
          else if (other2==other3)
            AddAtom(1,bond_ref,other1,other2,bond,120.0,180.0);
        } 
      } 
    }
  }

  // Obtain cartesian representation
  Z2C(0); 
}

