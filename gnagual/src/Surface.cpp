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

#include <GNagual.h>
#include <Surface.h>
#include <System.h>
#include <Grid.h>
#include <Spin.h>
#include <Set.h>
#include <Shell.h>
#include <Plotter.h>
#include <Evaluator.h>
#include <Math.h>

using namespace std;

Surface::Surface(System* s,Grid *g)
{
  sys = s;
  grid = g;
  scalar = new double[grid->npoint[0]];
};

void Surface::GetCubeValues(int i,int j,int k, double* val)
{
  int l,ll;
  ll = grid->npoint[1]*grid->npoint[2]*(k-1);

  l = ll + grid->npoint[1]*j + i-1;    
  val[0] = scalar[l];
  val[1] = scalar[l+1];
  l = ll + grid->npoint[1]*(j-1) + i;
  val[2] = scalar[l];
  val[3] = scalar[l-1];
  ll = grid->npoint[1]*grid->npoint[2]*k;
  l = ll + grid->npoint[1]*j + i-1;
  val[4] = scalar[l];
  val[5] = scalar[l+1];
  l = ll + grid->npoint[1]*(j-1) + i;
  val[6] = scalar[l];
  val[7] = scalar[l-1];
}

void Surface::Setup(void)
{

}

void Surface::Build(void)
{
  if (type==COLOURED_ISO)
  {
    BuildDensity();
  }
  else
  {
    if ( strncasecmp(name,"ORBITAL",7) == 0 )  
      BuildOrbital();
    else if ( strncasecmp(name,"DENSITY",7) == 0 ) 
      BuildDensity();
    else if ( strncasecmp(name,"SPECIES DENSITY",15) == 0 ) 
      BuildSpeciesDensity();
    else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
      BuildChargeDensity();
    else if ( strncasecmp(name,"SPIN DENSITY",12) == 0 ) 
      BuildElectronicSpinDensity();
    else if ( strncasecmp(name,"SHANNON ENTROPY",15) == 0 ) 
      BuildShannon();
    else if ( strncasecmp(name,"ESP",3) == 0 || strstr(name,"(ESP)") != NULL ) 
      BuildElectrostaticPotential();
    else if ( strncasecmp(name,"GTO",3) == 0 ) 
      BuildGTO();
  }
};

void Surface::BuildGTO()
{
  int i,ig,ib,is,io,j,k,ni,nj,nk;
  ig = set_number;
  if (ig>=sys->nspin) 
  {
    cout << "WARNING: Attempt to plot for non existing species"<<endl;
    return;
  } 
  else if (orbital_number>=sys->spin[ig]->nco)
  {
    cout << "WARNING: Attempt to non existing orbital"<<endl;
    return;
  }

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  nscalar = 0;
  i = 0;
  for (ib=0;ib<(signed)sys->spin[ig]->basis.size();ib++)
  {
    for (is=0;is<(signed)sys->spin[ig]->basis[ib]->shell.size();is++)
    {
      Shell *shell = sys->spin[ig]->basis[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==orbital_number) 
        {
          double x0,y0,z0;
          x0 = sys->spin[ig]->basis[ib]->x;
          y0 = sys->spin[ig]->basis[ib]->y;
          z0 = sys->spin[ig]->basis[ib]->z;
          for ( k = 0 ; k < nk ; k++ )
          {
            for ( j = 0 ; j < nj ; j++ )
            {
              grid->GetRowOfPoints(1,j,k,x,y,z);
              for (int ii=0;ii<ni;ii++) 
              {
                x[ii] -= x0;
                y[ii] -= y0;
                z[ii] -= z0;
              }
              shell->EvaluateGTO(x,y,z,&scalar[nscalar],io,ni);
              nscalar += ni;
            }
          }
          return;
        }
        i++;
      }
    }
  }
};

void Surface::BuildOrbital()
{
  int i,ig,io,j,k,ni,nj,nk;

  ig = set_number;
  if (ig>=sys->nspin) 
  {
    cout << "WARNING: Attempt to plot for non existing species"<<endl;
    return;
  } 
  else if (orbital_number>=sys->spin[ig]->nco)
  {
    cout << "WARNING: Attempt to non existing orbital"<<endl;
    return;
  }
  

  int nco = sys->spin[ig]->nco;
  double gto[nco];
  double moc[nco];
  sys->spin[ig]->C->GetRowValues( orbital_number , moc );

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  nscalar = 0;
  for ( k = 0 ; k < nk ; k++ )
  {
    for ( j = 0 ; j < nj ; j++ )
    {
      grid->GetRowOfPoints(1,j,k,x,y,z);
      for (i=0;i<ni;i++)
      {
        GetBasis(sys->spin[ig]->basis,x[i],y[i],z[i],gto); 
        scalar[nscalar] = 0.0;
        for (io=0;io<nco;io++)
          scalar[nscalar] += moc[io]*gto[io];
        nscalar++;
      }
    }
  }
};



void Surface::GetDerivatives(double* x, double *y, double *z, double *dx,
double *dy, double *dz, int n)
{
  if (type==COLOURED_ISO)
  {
    GetDensityDerivatives(sys,0,x,y,z,dx,dy,dz,n);
  }
  else
  {
    if ( strncasecmp(name,"ORBITAL",7) == 0 ) 
      GetOrbitalDerivatives(sys->spin[set_number],orbital_number,x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"DENSITY",7) == 0 ) 
      GetDensityDerivatives(sys,perturbation_number,x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"SPECIES DENSITY",15) == 0 ) 
      GetSpeciesDensityDerivatives(sys->spin[set_number],perturbation_number,
                                   x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
      GetChargeDensityDerivatives(sys,perturbation_number,x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"SPIN DENSITY",12) == 0 ) 
      GetElectronicSpinDensityDerivatives(sys,perturbation_number,x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"SHANNON ENTROPY",15) == 0 ) 
      GetShannonDerivatives(sys,perturbation_number,x,y,z,dx,dy,dz,n);
    else if ( strncasecmp(name,"ESP",3) == 0 || strstr(name,"(ESP)") != NULL ) 
      cout << "Derivatives NOT available for ESP " <<endl;
    else if ( strncasecmp(name,"GTO",3) == 0 ) 
      GetGTODerivatives(sys->spin[set_number]->basis,orbital_number,x,y,z,dx,dy,dz,n);
  }
};


void Surface::GetValues(double* x, double *y, double *z, double *v, int n)
{
  if ( strncasecmp(name,"ORBITAL",7) == 0 ) 
    GetOrbital(sys->spin[set_number],orbital_number,x,y,z,v,n);
  else if ( strncasecmp(name,"DENSITY",7) == 0 )  
    GetDensity(sys,perturbation_number,x,y,z,v,n);
  else if ( strncasecmp(name,"SPECIES DENSITY",15) == 0 )  
    GetSpeciesDensity(sys->spin[set_number],perturbation_number,x,y,z,v,n);
  else if ( strncasecmp(name,"CHARGE DENSITY",14) == 0 ) 
    GetChargeDensity(sys,perturbation_number,x,y,z,v,n);
  else if ( strncasecmp(name,"SPIN DENSITY",12) == 0 ) 
    GetElectronicSpinDensity(sys,perturbation_number,x,y,z,v,n);
  else if ( strncasecmp(name,"SHANNON ENTROPY",15) == 0 ) 
    GetShannon(sys,perturbation_number,x,y,z,v,n);
  else if ( strncasecmp(name,"GTO",3) == 0 ) 
    GetGTO(sys->spin[set_number]->basis,orbital_number,x,y,z,v,n);
};

void Surface::GetExtrema(double *vmin, double *vmax)
{
  int i;
  double val;

  *vmin = 99999999.0;
  *vmax = 0.0;
  for (i=0;i<nscalar;i++)
  {
    val = GN_ABS(scalar[i]);
    *vmin = GN_MIN(*vmin,val);
    *vmax = GN_MAX(*vmax,val);
  }
}

// Functions

// Density like functions
void Surface::BuildDensityLike(Matrix *P, vector<Set*> basis)
{
  int i,io,j,jo,k,ni,nj,nk;

  int nco = basis[basis.size()-1]->ul + 1;
  double gto[nco];
  double row[nco];

  ni = grid->npoint[1];
  nj = grid->npoint[2];
  nk = grid->npoint[3];
  double x[ni];
  double y[ni];
  double z[ni];

  if (!P) 
  {
    cout << "Required matrix is not available" << endl;
    return;
  }

  nscalar = 0;
  for ( k = 0 ; k < nk ; k++ )
  {
    for ( j = 0 ; j < nj ; j++ )
    {
      grid->GetRowOfPoints(1,j,k,x,y,z);
      for (i=0;i<ni;i++)
      {
        GetBasis(basis,x[i],y[i],z[i],gto); 
        scalar[nscalar] = 0.0;
        for (io=0;io<nco;io++)
        {
          P->GetRowValues( io , row );
          for (jo=0;jo<nco;jo++)
            scalar[nscalar] += row[jo]*gto[io]*gto[jo];
        }
        nscalar++;
      }
    }
  }
};

// Density
void Surface::BuildDensity()
{
  int i,ig,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  for (ig=0;ig<sys->nspin;ig++)
  {
    set_number = ig;
    BuildSpeciesDensity();
    for (i=0;i<nbuf;i++)
      buf[i] += scalar[i];
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i];
}

// Species density
void Surface::BuildSpeciesDensity()
{
  if (set_number>=sys->nspin) 
  {
    cout << "WARNING: Attempt to plot for non existing species"<<endl;
    return;
  } 
  BuildDensityLike(sys->spin[set_number]->P[perturbation_number],
                   sys->spin[set_number]->basis);
};

// Charge density
void Surface::BuildChargeDensity()
{
  int i,ig,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  for (ig=0;ig<sys->nspin;ig++)
  {
    set_number = ig;
    BuildSpeciesDensity();
    for (i=0;i<nbuf;i++)
      buf[i] += scalar[i]*sys->spin[ig]->charge;
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i];
}

// Electronic density
void Surface::BuildElectronicSpinDensity()
{
  int i,ig,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  bool leader = true;
  for (ig=0;ig<sys->nspin;ig++)
  {
    if ( strncasecmp(sys->spin[ig]->name,"ELECT",5) == 0 ) 
    {
      set_number = ig;
      BuildSpeciesDensity();
      if (leader)
      {
        for (i=0;i<nbuf;i++)
          buf[i] += scalar[i];
        leader = false;
      }
      else
      {
        for (i=0;i<nbuf;i++)
          buf[i] -= scalar[i];
      }
   }
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i];
}

void Surface::BuildShape()
{
  int i,ig,nbuf;
  
  nbuf = grid->npoint[0];
  double buf[nbuf];

  for (i=0;i<nbuf;i++)
    buf[i] = 0.0;

  int nelec = 0;
  for (ig=0;ig<sys->nspin;ig++)
  {
    if ( strncasecmp(sys->spin[ig]->name,"ELECT",5) == 0 ) 
    {
      nelec += sys->spin[ig]->npart;
      set_number = ig;
      BuildSpeciesDensity();
      for (i=0;i<nbuf;i++)
        buf[i] += scalar[i];
    }
  }

  for (i=0;i<nbuf;i++)
    scalar[i] = buf[i]/double(nelec);
}

// Shannon
void Surface::BuildShannon()
{
  int i;
  double shape;

  BuildShape();

  for (i=0;i<nscalar;i++)
  {
    shape = scalar[i];
    scalar[i] = -shape*log(shape);
  }
}

// ===================================
// Electrostatic potential
void Surface::BuildElectrostaticPotential()
{
}

