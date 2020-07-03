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

#include <Evaluator.h>
#include <System.h>
#include <Molecule.h>
#include <Atom.h>
#include <Spin.h>
#include <Set.h>
#include <Shell.h>
#include <Math.h>

using namespace std;

void GetGTODerivatives(vector<Set*> basis, int id, double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,ib,is,io;

  double wx[n];
  double wy[n];
  double wz[n];

  i = 0;
  for (ib=0;ib<(signed)basis.size();ib++)
  {
    for (is=0;is<(signed)basis[ib]->shell.size();is++)
    {
      Shell *shell = basis[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==id) 
        {
          double x0,y0,z0;

          x0 = basis[ib]->x;
          y0 = basis[ib]->y;
          z0 = basis[ib]->z;
          for (int ii=0;ii<n;ii++) 
          {
            wx[ii] = x[ii] - x0;
            wy[ii] = y[ii] - y0;
            wz[ii] = z[ii] - z0;
          }
          shell->EvaluateGTODerivative(wx,wy,wz,dx,dy,dz,io,n);
          return;
        }
        i++;
      }
    }
  }
};

void GetBasisDerivatives(vector<Set*> basis, double x, double y, double z,
 double *dx, double *dy, double *dz)
{
  int i,ib,is,io,lla;
  Shell *shell;

  double xx[1];
  double yy[1];
  double zz[1];
  double dxx[1];
  double dyy[1];
  double dzz[1];

  for (ib=0;ib<(signed)basis.size();ib++)
  {
    xx[0] = x - basis[ib]->x;
    yy[0] = y - basis[ib]->y;
    zz[0] = z - basis[ib]->z;
    for (is=0;is<(signed)basis[ib]->shell.size();is++)
    {
      shell = basis[ib]->shell[is];
      lla = shell->ll + basis[ib]->ll;
      for (io=0;io<shell->nco;io++)
      {
        i = lla + io;
        shell->EvaluateGTODerivative(xx,yy,zz,dxx,dyy,dzz,io,1);
        dx[i] = dxx[0];
        dy[i] = dyy[0];
        dz[i] = dzz[0];
      }
    }
  }
};


void GetOrbitalDerivatives(Spin* spin, int id, double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,io;

  int nco = spin->nco;
  double moc[nco];
  spin->C->GetRowValues( id , moc );

  double gto_dx[nco];
  double gto_dy[nco];
  double gto_dz[nco];

  for (i=0;i<n;i++) 
  {
    GetBasisDerivatives(spin->basis,x[i],y[i],z[i],gto_dx,gto_dy,gto_dz);
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
    for (io=0;io<nco;io++) 
    {
      dx[i] += gto_dx[io]*moc[io];
      dy[i] += gto_dy[io]*moc[io];
      dz[i] += gto_dz[io]*moc[io];
    }
  }
};



void GetDensityDerivatives(System* sys,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,ig;

  double dxbuf[n],dybuf[n],dzbuf[n];

  for (i=0;i<n;i++) 
  {
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
  }
  for (ig=0;ig<sys->nspin;ig++)
  {
    GetSpeciesDensityDerivatives(sys->spin[ig],pid,x,y,z,dxbuf,dybuf,dzbuf,n);
    for (i=0;i<n;i++) 
    {
      dx[i] += dxbuf[i];
      dy[i] += dybuf[i];
      dz[i] += dzbuf[i];
    }
  }
}

void GetChargeDensityDerivatives(System* sys,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,ig;

  double dxbuf[n],dybuf[n],dzbuf[n];

  for (i=0;i<n;i++) 
  {
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
  }
  for (ig=0;ig<sys->nspin;ig++)
  {
    GetSpeciesDensityDerivatives(sys->spin[ig],pid,x,y,z,dxbuf,dybuf,dzbuf,n);
    for (i=0;i<n;i++) 
    {
      dx[i] += dxbuf[i]*sys->spin[ig]->charge;
      dy[i] += dybuf[i]*sys->spin[ig]->charge;
      dz[i] += dzbuf[i]*sys->spin[ig]->charge;
    }
  }
}

void GetElectronicSpinDensityDerivatives(System* sys,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
}

void GetShapeDerivatives(System* sys,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i,ig,nelec;

  double dxbuf[n],dybuf[n],dzbuf[n];

  for (i=0;i<n;i++) 
  {
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
  }
  nelec = 0;
  for (ig=0;ig<sys->nspin;ig++)
  {
    if ( strncasecmp(sys->spin[ig]->name,"ELECT",5) == 0 ) 
    {
      nelec += sys->spin[ig]->npart;
      GetSpeciesDensityDerivatives(sys->spin[ig],pid,x,y,z,dxbuf,dybuf,dzbuf,n);
      for (i=0;i<n;i++) 
      {
        dx[i] += dxbuf[i];
        dy[i] += dybuf[i];
        dz[i] += dzbuf[i];
      }
    }
  }
  if ( nelec == 0 ) return;
  for (i=0;i<n;i++) 
  {
    dx[i] /= double(nelec);
    dy[i] /= double(nelec);
    dz[i] /= double(nelec);
  }
}
      
void GetGTO(vector<Set*> basis, int id, double* x, double *y, double *z, double *v, int n)
{
  int i,ib,is,io;

  double wx[n];
  double wy[n];
  double wz[n];

  i = 0;
  for (ib=0;ib<(signed)basis.size();ib++)
  {
    for (is=0;is<(signed)basis[ib]->shell.size();is++)
    {
      Shell *shell = basis[ib]->shell[is];
      for (io=0;io<shell->nco;io++)
      {
        if (i==id) 
        {
          double x0,y0,z0;

          x0 = basis[ib]->x;
          y0 = basis[ib]->y;
          z0 = basis[ib]->z;
          for (int ii=0;ii<n;ii++) 
          {
            wx[ii] = x[ii] - x0;
            wy[ii] = y[ii] - y0;
            wz[ii] = z[ii] - z0;
          }
          shell->EvaluateGTO(wx,wy,wz,v,io,n);
          return;
        }
        i++;
      }
    }
  }
};

void GetBasis(vector<Set*> basis, double x, double y, double z,double *v)
{
  int i,ib,is,io,lla;

  double xx[1];
  double yy[1];
  double zz[1];

  for (ib=0;ib<(signed)basis.size();ib++)
  {
    for (is=0;is<(signed)basis[ib]->shell.size();is++)
    {
      Shell *shell = basis[ib]->shell[is];
      lla = shell->ll + basis[ib]->ll;
      for (io=0;io<shell->nco;io++)
      {
        i = lla + io;
        if (basis[ib]->is_neighbor)
        {
          xx[0] = x - basis[ib]->x;
          yy[0] = y - basis[ib]->y;
          zz[0] = z - basis[ib]->z;
          shell->EvaluateGTO(xx,yy,zz,&v[i],io,1);
        }
        else v[i] = 0.0;
      }
    }
  }
};


void GetOrbital(Spin* spin, int id, double* x, double *y, double *z, double *v, int n)
{
  int i,io;

  int nco = spin->nco;
  double moc[nco];
  spin->C->GetRowValues( id , moc );

  double gto[nco];

  for (i=0;i<n;i++) 
  {
    GetBasis(spin->basis,x[i],y[i],z[i],gto);
    v[i] = 0.0;
    for (io=0;io<nco;io++) 
      v[i] += gto[io]*moc[io];
  }
};


void GetDensity(System* sys,int pid,double* x, double *y, double *z, double *v, int n)
{
  int i,ig;
  
  double buf[n];

  for (i=0;i<n;i++)
    v[i] = 0.0;
  for (ig=0;ig<sys->nspin;ig++)
  {
    GetSpeciesDensity(sys->spin[ig],pid,x,y,z,buf,n);
    for (i=0;i<n;i++)
      v[i] += buf[i];
  }
}


void GetChargeDensity(System* sys,int pid,double* x, double *y, double *z, double *v, int n)
{
  int i,ig;
  
  double buf[n];

  for (i=0;i<n;i++)
    v[i] = 0.0;
  for (ig=0;ig<sys->nspin;ig++)
  {
    GetSpeciesDensity(sys->spin[ig],pid,x,y,z,buf,n);
    for (i=0;i<n;i++)
      v[i] += buf[i]*sys->spin[ig]->charge;
  }
}


void GetElectronicSpinDensity(System* sys,int pid, double* x, double *y, double *z, double *v, int n)
{
  int i,ig;
  
  double buf[n];

  bool leader = true;
  for (ig=0;ig<sys->nspin;ig++)
  {
    if ( strncasecmp(sys->spin[ig]->name,"ELECT",5) == 0 ) 
    {
      GetSpeciesDensity(sys->spin[ig],pid,x,y,z,buf,n);
      if (leader)
      {
        for (i=0;i<n;i++)
          v[i] = buf[i];
        leader = false;
      }
      else
      {
        for (i=0;i<n;i++)
          v[i] -= buf[i];
      }
    }
  }
}

void GetShape(System* sys,int pid,double* x, double *y, double *z, double *v, int n)
{
  int i,ig,nelec;
  
  double buf[n];

  for (i=0;i<n;i++)
    v[i] = 0.0;
  nelec = 0;
  for (ig=0;ig<sys->nspin;ig++)
  {
    if ( strncasecmp(sys->spin[ig]->name,"ELECT",5) == 0 ) 
    {
      nelec += sys->spin[ig]->npart;
      GetSpeciesDensity(sys->spin[ig],pid,x,y,z,buf,n);
      for (i=0;i<n;i++)
        v[i] += buf[i];
    }
  }
  if (nelec==0) return;
  for (i=0;i<n;i++)
    v[i] /= double(nelec);
}

// Functions

// ==========================
// Density like functions
void GetDensityLike(Matrix *P, vector<Set*> basis,
double* x, double *y, double *z, double *v, int n)
{
  int i,io,jo;

  int nco = basis[basis.size()-1]->ul + 1;
  double s;
  double row[nco];
  double gto[nco];

  if (!P) 
  {
    cout << "Required matrix is not available" << endl;
    return;
  }

  for (int i=0;i<n;i++) 
  {
    GetBasis(basis,x[i],y[i],z[i],gto); 
    v[i] = 0.0;
    for (io=0;io<nco;io++) 
    {
      P->GetRowValues( io , row );
      s = row[io]*gto[io];
      for (jo=0;jo<io;jo++) 
        s += 2.0*row[jo]*gto[jo];
      v[i] += s*gto[io];
    }
  }
};

void GetDensityLikeDerivatives(Matrix *P,vector<Set*> basis,
double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int io,jo;

  int nco = basis[basis.size()-1]->ul + 1;
  double row[nco];

  double gto[nco];
  double gto_dx[nco];
  double gto_dy[nco];
  double gto_dz[nco];

  if (!P) 
  {
    cout << "Required matrix is not available" << endl;
    return;
  }

  for (int i=0;i<n;i++) 
  {
    GetBasis(basis,x[i],y[i],z[i],gto); 
    GetBasisDerivatives(basis,x[i],y[i],z[i],gto_dx,gto_dy,gto_dz);
    dx[i] = 0.0;
    dy[i] = 0.0;
    dz[i] = 0.0;
    for (io=0;io<nco;io++) 
    {
      P->GetRowValues( io , row );
      for (jo=0;jo<nco;jo++) 
      {
        dx[i] += row[jo]*(gto_dx[io]*gto[jo]+gto[io]*gto_dx[jo]);
        dy[i] += row[jo]*(gto_dy[io]*gto[jo]+gto[io]*gto_dy[jo]);
        dz[i] += row[jo]*(gto_dz[io]*gto[jo]+gto[io]*gto_dz[jo]);
      }
    }
  }
};

// ==========================
// Species density
void GetSpeciesDensity(Spin* spin,int pid,double* x, double *y, double *z, double *v, int n)
{
  GetDensityLike(spin->P[pid],spin->basis,x,y,z,v,n);
};

void GetSpeciesDensityDerivatives(Spin* spin,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  GetDensityLikeDerivatives(spin->P[pid],spin->basis,x,y,z,dx,dy,dz,n);
};

// ==========================
// Shannon entropy based on shape function
void GetShannon(System* sys, int pid,double* x, double *y, double *z, double *v, int n)
{
  int i;
  double shape;

  GetShape(sys,pid,x,y,z,v,n);

  for (i=0;i<n;i++)
  {
    shape = v[i];
    v[i] = -shape*log(shape);
  }
}

void GetShannonDerivatives(System* sys,int pid,double* x, double *y, double *z, double *dx, double *dy, double *dz, int n)
{
  int i;
  double shape,dshape;
  double v[n];

  GetShape(sys,pid,x,y,z,v,n);
  GetShapeDerivatives(sys,pid,x,y,z,dx,dy,dz,n);

  for (i=0;i<n;i++) 
  {
    shape = v[i];
    dshape = dx[i]; dx[i] = - log(shape) - dshape;
    dshape = dy[i]; dy[i] = - log(shape) - dshape;
    dshape = dz[i]; dz[i] = - log(shape) - dshape;
  }
}
      
