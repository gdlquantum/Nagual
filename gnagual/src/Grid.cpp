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
#include <vector>

#include <Parameter.h>
#include <Grid.h>
#include <Vector.h>

using namespace std;

Grid::Grid()
{
  double lat[12] = {-GN_GS,-GN_GS,-GN_GS,
                     GN_GS,-GN_GS,-GN_GS,
                    -GN_GS, GN_GS,-GN_GS, 
                    -GN_GS,-GN_GS, GN_GS};
  SetLattice(lat,0.02*GN_GS);
}



void Grid::GetCubePoints(int i,int j,int k, Vector r[8])
{
  int l,ll;
  ll = npoint[1]*npoint[2]*(k-1);

  r[0] = Vector(vertex[0]+dstep[1]*(i-1)+dstep[2]*j    +dstep[3]*(k-1));
  r[1] = Vector(vertex[0]+dstep[1]*i    +dstep[2]*j    +dstep[3]*(k-1));
  r[2] = Vector(vertex[0]+dstep[1]*i    +dstep[2]*(j-1)+dstep[3]*(k-1));
  r[3] = Vector(vertex[0]+dstep[1]*(i-1)+dstep[2]*(j-1)+dstep[3]*(k-1));
  r[4] = Vector(vertex[0]+dstep[1]*(i-1)+dstep[2]*j    +dstep[3]*k);
  r[5] = Vector(vertex[0]+dstep[1]*i    +dstep[2]*j    +dstep[3]*k);
  r[6] = Vector(vertex[0]+dstep[1]*i    +dstep[2]*(j-1)+dstep[3]*k);
  r[7] = Vector(vertex[0]+dstep[1]*(i-1)+dstep[2]*(j-1)+dstep[3]*k);
}


void Grid::GetLattice(double* lat)
{
  lat[0] = vertex[0].x;
  lat[1] = vertex[0].y;
  lat[2] = vertex[0].z;
  for (int i=1;i<4;i++)
  {
    lat[3*i] = vertex[i].x;
    lat[3*i+1] = vertex[i].y;
    lat[3*i+2] = vertex[i].z;
  }
}


void Grid::SetLattice(double* lat, double boxmesh)
{
  for (int i=0;i<4;i++)
  {
    vertex[i] = Vector(&lat[3*i]);
    if (i>0) 
    {
      dstep[i] = vertex[i]-vertex[0];
      double norm = dstep[i].Norm();
      npoint[i] = int(norm/boxmesh) + 1;
      dstep[i] ^= boxmesh;
    }
  }
  npoint[0] = npoint[1]*npoint[2]*npoint[3];
}

void Grid::GetRowOfPoints(int indx, int b,int c, double* x,double *y, double *z)
{
  int i,n;

  n = npoint[indx];

  Vector v = vertex[0];

  if (indx==1)
  {
    v += dstep[2]*b;
    v += dstep[3]*c;
    for (i=0;i<n;i++)
    {
      x[i] = v.x;
      y[i] = v.y;
      z[i] = v.z; 
      v += dstep[1];
    }
  }
  else if (indx==2)
  {
    v += dstep[3]*b;
    v += dstep[1]*c;
    for (i=0;i<n;i++)
    {
      x[i] = v.x;
      y[i] = v.y;
      z[i] = v.z; 
      v += dstep[2];
    }
  }
  else if (indx==3)
  {
    v += dstep[1]*b;
    v += dstep[2]*c;
    for (i=0;i<n;i++)
    {
      x[i] = v.x;
      y[i] = v.y;
      z[i] = v.z; 
      v += dstep[3];
    }
  }
}
