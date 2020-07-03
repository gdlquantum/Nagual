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
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <Shell.h>
#include <Math.h>
#include <Vector.h>

using namespace std;

Shell::Shell(void)
{
  z.clear();
  d.clear();

  is_normalized = false;
  radius = 1.0/GN_TOL_NUM;
};

void Shell::Normalize(void)
{
  if (is_normalized) return;

  int i,j,n;
  double factor,norm;

  n = z.size();

  // Normalized GTOs
  norm = pow(2.0,double(l))*pow(2.0/GN_PI,0.75);
  for ( i = 0 ; i < n ; i++ )
  {
    d[i] = norm*pow(z[i],(2.0*l+3.0)*0.25)*d[i];
  }

  factor = pow(GN_PI,1.5)/pow(2.0,(double)l);
  norm = 0.0;
  for ( i = 0 ; i < n ; i++ )
  {
    for ( j = 0 ; j < n ; j++ )
    {
      norm += d[i]*d[j]*factor/pow(z[i]+z[j],l+1.5);
    }
  }

  // Normalized STOs
  int lx,ly,lz;
  nco = 0;
  for ( lx = l ; lx >= 0 ; lx-- )
  {
    for ( ly = l - lx ; ly >= 0 ; ly-- )
    {
      lz = l - lx - ly;
      ncsto.push_back( 1.0/sqrt(norm*DoubleFactorial(2*lx-1)*
                   DoubleFactorial(2*ly-1)*DoubleFactorial(2*lz-1)) );
      nco++;
    }
  }

  is_normalized = true;
}

void Shell::Print(char* filename)
{
  ofstream f(filename,ios::app);

  double norm = pow(2.0,double(l))*pow(2.0/GN_PI,0.75);
  int k = z.size();
  f << setw(2) << fixed << n << " "
    << setw(2) << fixed << l << " "
    << setw(3) << fixed << k <<endl;
  for ( int i = 0 ; i < k ; i++ )
  {
     f << setw(20) << fixed << setprecision(10) << z[i] << "     ";
     if (type==PRIMARY_BASIS_SHELL)
       f << setw(20) << fixed << setprecision(10)<< d[i]/(norm*pow(z[i],(2.0*l+3.0)*0.25))<<endl;
     else f << setw(20) << fixed << setprecision(10)<< d[i] <<endl;
  }
  f.close();
}


void Shell::EvaluateGTO(double *xc, double *yc, double *zc, double* res,int indx, int nr)
{
  Vector a;
  int i,ig,lx,ly,lz,n;
  double expf,r2;
  double xyzp;

  n = z.size();
  for ( ig = 0 ; ig < nr ; ig++ )
  {
    a = Vector(xc[ig],yc[ig],zc[ig]);
    r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (sqrt(r2)<radius)
    {
      expf = 0.0;
      for ( i = 0 ; i < n ; i++ )
      {
        expf += d[i]*exp(-z[i]*r2);
      }
      int ico = 0;
      for ( lx = l ; lx >= 0 ; lx-- )
      {
        for ( ly = l - lx ; ly >= 0 ; ly-- )
        {
          lz = l - lx - ly;
          if (ico==indx)
          {
            xyzp = 1.0;
            for (i=1;i<=lx;i++) xyzp *= a.x;
            for (i=1;i<=ly;i++) xyzp *= a.y;
            for (i=1;i<=lz;i++) xyzp *= a.z;
            res[ig] = ncsto[ico]*xyzp*expf;
          }
          ico++;
        }
      }
    } else res[ig] = 0.0;
  }
};


void Shell::EvaluateGTODerivative(double *xc, double *yc, double *zc, double* dx, double* dy, double *dz,int indx, int nr)
{
  Vector a;
  int i,ig,lx,ly,lz,n;
  double expf,expfp,r2,t;
  double rxp[GN_MAX_L+2],ryp[GN_MAX_L+2],rzp[GN_MAX_L+2];

  n = z.size();
  for ( ig = 0 ; ig < nr ; ig++ )
  {
    a = Vector(xc[ig],yc[ig],zc[ig]);
    r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    expf = 0.0;
    expfp = 0.0;
    for ( i = 0 ; i < n ; i++ )
    {
      t = d[i]*exp(-z[i]*r2);
      expf += t;
      expfp += t*2.0*z[i];
    }
    rxp[0] = 1.0;
    for (lx=1;lx<=l+1;lx++) rxp[lx] = rxp[lx-1]*a.x;
    ryp[0] = 1.0;
    for (ly=1;ly<=l+1;ly++) ryp[ly] = ryp[ly-1]*a.y;
    rzp[0] = 1.0;
    for (lz=1;lz<=l+1;lz++) rzp[lz] = rzp[lz-1]*a.z;
    int ico = 0;
    for ( lx = l ; lx >= 0 ; lx-- )
    {
      for ( ly = l - lx ; ly >= 0 ; ly-- )
      {
        lz = l - lx - ly;
        if (ico==indx)
        {
          t = ncsto[ico]*ryp[ly]*rzp[lz];
          dx[ig] = t*expfp*rxp[lx+1];
          if (lx>0) dx[ig] -= double(lx)*t*expf*rxp[lx-1];
          t = ncsto[ico]*rxp[lx]*rzp[lz];
          dy[ig] = t*expfp*ryp[ly+1];
          if (ly>0) dy[ig] -= double(ly)*t*expf*ryp[ly-1];
          t = ncsto[ico]*ryp[ly]*rxp[lx];
          dz[ig] = t*expfp*rzp[lz+1];
          if (lz>0) dz[ig] -= double(lz)*t*expf*rzp[lz-1];
        }
        ico++;
      }
    }
  }
};

void Shell::EvaluateGTOHessian(double *xc, double *yc, double *zc, double* dxx, double* dxy, double *dxz,double *dyy, double *dyz, double *dzz, int indx, int nr)
{
  Vector a;
  int i,ig,lx,ly,lz,n;
  double expf,expfp,expfpp,r2,t;
  double rxp[GN_MAX_L+3],ryp[GN_MAX_L+3],rzp[GN_MAX_L+3];

  n = z.size();
  for ( ig = 0 ; ig < nr ; ig++ )
  {
    a = Vector(xc[ig],yc[ig],zc[ig]);
    r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    expf = 0.0;
    expfp = 0.0;
    expfpp = 0.0;
    for ( i = 0 ; i < n ; i++ )
    {
      t = d[i]*exp(-z[i]*r2);
      expf += t;
      expfp += t*2.0*z[i];
      expfpp += t*4.0*z[i]*z[i];
    }
    rxp[0] = 1.0;
    for (lx=1;lx<=l+2;lx++) rxp[lx] = rxp[lx-1]*a.x;
    ryp[0] = 1.0;
    for (ly=1;ly<=l+2;ly++) ryp[ly] = ryp[ly-1]*a.y;
    rzp[0] = 1.0;
    for (lz=1;lz<=l+2;lz++) rzp[lz] = rzp[lz-1]*a.z;
    int ico = 0;
    for ( lx = l ; lx >= 0 ; lx-- )
    {
      for ( ly = l - lx ; ly >= 0 ; ly-- )
      {
        lz = l - lx - ly;
        if (ico==indx)
        {
          // xx
          t = ncsto[ico]*ryp[ly]*rzp[lz];
          dxx[ig] = t*expfpp*rxp[lx+2];
          dxx[ig] -= (1,0+2.0*double(lx))*t*expfp*rxp[lx];
          if (lx>1) dxx[ig] += double(lx*(lx-1))*t*expf*rxp[lx-2];
          // xy
          t = ncsto[ico]*rzp[lz];
          dxy[ig] = t*expfpp*rxp[lx+1]*ryp[ly+1];
          if (lx>0) dxy[ig] -= double(lx)*t*expfp*rxp[lx-1]*ryp[ly];
          if (ly>0) dxy[ig] -= double(ly)*t*expfp*rxp[lx]*ryp[ly-1];
          if (lx>0&&ly>0) dxy[ig] += double(lx*ly)*t*expf*rxp[lx-1]*ryp[ly-1];
          // xz
          t = ncsto[ico]*ryp[ly];
          dxz[ig] = t*expfpp*rxp[lx+1]*rzp[lz+1];
          if (lx>0) dxz[ig] -= double(lx)*t*expfp*rxp[lx-1]*rzp[lz];
          if (lz>0) dxz[ig] -= double(lz)*t*expfp*rxp[lx]*rzp[lz-1];
          if (lx>0&&lz>0) dxz[ig] += double(lx*lz)*t*expf*rxp[lx-1]*rzp[lz-1];
          // yy
          t = ncsto[ico]*rxp[lx]*rzp[lz];
          dyy[ig] = t*expfpp*ryp[ly+2];
          dyy[ig] -= (1,0+2.0*double(ly))*t*expfp*ryp[ly];
          if (ly>1) dyy[ig] += double(ly*(ly-1))*t*expf*ryp[ly-2];
          // yz
          t = ncsto[ico]*rxp[lx];
          dyz[ig] = t*expfpp*ryp[ly+1]*rzp[lz+1];
          if (ly>0) dyz[ig] -= double(ly)*t*expfp*ryp[ly-1]*rzp[lz];
          if (lz>0) dyz[ig] -= double(lz)*t*expfp*ryp[ly]*rzp[lz-1];
          if (ly>0&&lz>0) dyz[ig] += double(ly*lz)*t*expf*ryp[ly-1]*rzp[lz-1];
          // zz
          t = ncsto[ico]*rxp[lx]*ryp[ly];
          dzz[ig] = t*expfpp*rzp[lz+2];
          dzz[ig] -= (1,0+2.0*double(lz))*t*expfp*rzp[lz];
          if (lz>1) dzz[ig] += double(lz*(lz-1))*t*expf*rzp[lz-2];
        }
        ico++;
      }
    }
  }
};

// Compute Radial potential (only exponential part)
double Shell::RadialPotential(double r)
{
  int i;
  double expf,r2;

  // Build potential
  r2 = r*r;

  // Contract the potential of the current shell 
  expf = 0.0;
  for (i=0;i<(signed)z.size();i++)
    expf += d[i]*exp(-z[i]*r2);

  return expf;
}

double Shell::Radius()
{
  return radius;
}

void Shell::EvaluateRadius()
{
  size_t i,imin;
  double delta,dg,g,ldt,zr;

  // Get smallest exponent
  imin = 0;
  for (i=1;i<z.size();i++)
  {
    if (z[i]<z[imin])
      imin = i;
  }

  ldt = log(GN_ABS(d[imin])/GN_ABS(1.0e-10));
  radius = GN_MAX(ldt/z[imin],GN_TOL_NUM);

  // For large radial powers the guess should be corrected ***
  if (l>0)
  {
    g = GN_MAX(sqrt(double(l)/(2.0*GN_ABS(d[imin]))),GN_TOL_NUM);
    if (g>radius) radius = 0.5*(radius+g);
  }

  for (i=0;i<100;i++)
  {
    zr = z[imin]*radius;
    g = ldt + double(l)*log(radius) - zr*radius;
    dg = double(l)/radius - 2.0*zr;
    delta = g/dg;
    radius = GN_MAX(radius-delta,GN_TOL_NUM);
    if (GN_ABS(delta)<1.0E-10) break;
  }
  if (GN_ABS(delta)>=1.0E-10) 
  {
    radius = 1.0/GN_TOL_NUM;
    cout << "Unable to find shell radius, infinity assumed" << endl;
  }
}


