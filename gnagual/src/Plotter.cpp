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
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <Parameter.h>
#include <Plotter.h>
#include <Viewer.h>
#include <DwPanel.h>
#include <GNagual.h>
#include <Vector.h>
#include <Grid.h>
#include <Surface.h>
#include <Math.h>
#include <Units.h>

char* stylenames[NPLOT_STYLES] = {"Solid","Points", "Translucent", "Mesh"};
char* typenames[NPLOT_TYPES] = {"Isosurface","Contours",
"Coloured Iso-Density","Coloured Plane","Curve"}; 

using namespace std;

Plotter::Plotter(GNagual *x)
{
  gn = x;
  surface.clear();
  grid = new Grid();
  type = ISO;
  style = SOLID_STYLE;
  use_lists = true;
  point_size = 3;
  plane_axis = 2;
  plane_point = 0;
}

int Plotter::NSurf()
{
  return (int) surface.size();
}

void Plotter::PlotISO( int n )
{
  if ( surface[n]->style == SOLID_STYLE )
  {
    glBegin( GL_TRIANGLES );
      Build(n);
    glEnd();
  }
  else if ( surface[n]->style == TRANSLUCENT_STYLE )
  {
    glEnable( GL_CULL_FACE );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    glBegin( GL_TRIANGLES );
      Build(n);
    glEnd();
    glDisable(GL_BLEND);
    glDisable( GL_CULL_FACE );
  }
  else if ( surface[n]->style == MESH_STYLE )
  {
    glLineWidth( surface[n]->line_width );
    glBegin( GL_LINES );
      Build(n);
    glEnd();
  }
  else if ( surface[n]->style == POINTS_STYLE )
  {
    glPointSize( surface[n]->point_size );
    glBegin( GL_POINTS );
      Build(n);
    glEnd();
  };
}

void Plotter::PlotCURVE( int n )
{
  glLineWidth( surface[n]->line_width );
  glDisable( GL_LIGHTING );
  glBegin( GL_LINES );
    Build(n); 
  glEnd();
  glEnable( GL_LIGHTING );
}

void Plotter::PlotCONTOUR( int n )
{
  double ccolor[3];

  glDisable( GL_LIGHTING );
  glLineWidth( surface[n]->line_width );
  surface[n]->iso = fmin;
  while ( surface[n]->iso < fmax )
  {
    RealToRGB( surface[n]->iso, ccolor );
    glColor3dv( ccolor );
    glBegin( GL_LINES );
    Build(n); 
    glEnd();
    surface[n]->iso += fstep;
  };
  glEnable( GL_LIGHTING );
}

void Plotter::PlotCOLOURED_ISO( int n )
{
  surface[n]->style = SOLID_STYLE;
  glEnable( GL_LIGHTING );
  glBegin( GL_TRIANGLES );
  Build(n); 
  glEnd();
}

void Plotter::PlotCOLOURED_PLANE( int n )
{
  surface[n]->style = SOLID_STYLE;
  glDisable( GL_LIGHTING );
  glBegin( GL_TRIANGLES );
  Build(n); 
  glEnd();
  glEnable( GL_LIGHTING );
}

void Plotter::GetColor( int n , double *color )
{
  if ( n < NSurf() )
  {
    for ( int i = 0 ; i < 4 ; i++ )
    {
      color[i] = surface[n]->color[i];
    };
  };
}

void Plotter::SetColor( int n , double *color )
{
  if ( n < NSurf() )
  {
    for ( int i = 0 ; i < 4 ; i++ )
    {
      surface[n]->color[i] = color[i];
    };
  };
}

void Plotter::RequestPlot(const QString& st,int sid,int oid, int pid,
double isovalue, bool setup)
{
  Surface pe(gn->sys,grid);
  strcpy( pe.name , (st.toLatin1()).data());
  pe.set_number = sid;
  pe.orbital_number = oid;
  pe.perturbation_number = pid;
  if ( strncmp(pe.name,"BOX",3) == 0 ) pe.type = CURVE; 
  else pe.type = type; 
  if ( pe.type == COLOURED_ISO ) pe.style = SOLID_STYLE; 
  else pe.style = style; 
  pe.shiny = shiny;
  pe.line_width = line_width;
  pe.point_size = point_size;

  if ( isovalue > 0 )
  {
    pe.color[0] = 1.0;
    pe.color[1] = 0.0;
    pe.color[2] = 0.0;
    pe.color[3] = 0.0;
  }
  else
  {
    pe.color[0] = 0.0;
    pe.color[1] = 0.0;
    pe.color[2] = 1.0;
    pe.color[3] = 0.0;
  }

  pe.iso = isovalue;  
  pe.list = 0; // No list available yet
  surface.push_back( new Surface( pe ) );
  int n = surface.size()-1;
  if (setup) surface[n]->Setup();
  else surface[n]->Build();
}

void Plotter::DoSurface( int n )
{
  if (n < NSurf()) 
  {
    GLfloat noneDir[4] = {0.2, 0.2, 0.2, 0.0};
    GLfloat whiteDir[4] = {2.0, 2.0, 2.0, 1.0};
    glColor3dv( surface[n]->color );
    if ( ! surface[n]->shiny )
    {
      glMaterialfv(GL_FRONT, GL_SPECULAR, noneDir);
    };
    if ( (! surface[n]->list ) || ( ! use_lists ) )
    {
      if (use_lists)
      {
        surface[n]->list = glGenLists( 1 );
        glNewList( surface[n]->list , GL_COMPILE );
      }
      if ( surface[n]->type == ISO ) PlotISO( n );
      else if ( surface[n]->type == CURVE ) PlotCURVE( n );
      else if ( surface[n]->type == CONTOUR ) PlotCONTOUR( n );
      else if ( surface[n]->type == COLOURED_ISO ) PlotCOLOURED_ISO( n );
      else if ( surface[n]->type == COLOURED_PLANE ) PlotCOLOURED_PLANE( n );
      if (use_lists) glEndList();
    };
    if ( use_lists ) glCallList( surface[n]->list );
    if ( ! surface[n]->shiny )
      glMaterialfv(GL_FRONT, GL_SPECULAR, whiteDir);
  }
}

void Plotter::Build( int n )
{
  int i,ip,it;
  int j,jgp,k,nt;
  Vector normal,pos;
  double rgb[3];
  double t[3][5][3],tn[3][5][3];

  // Active box
  int mini,minj,mink,maxi,maxj,maxk;
  mini = minj = mink = 0;
  maxi = grid->npoint[1];
  maxj = grid->npoint[2];
  maxk = grid->npoint[3];

  // Plotting
  if ( strncmp(surface[n]->name,"BOX",3) == 0 ) 
  {
    double LATTICE[12];

    grid->GetLattice(LATTICE);

    glEnable(GL_LIGHTING);
    glColor4d( 0.5 , 0.5 , 0.5 , 0.0 );
    glLineWidth( line_width );
    glEnd();
    glBegin( GL_LINE_STRIP );
    glVertex3d( LATTICE[0] , LATTICE[1] , LATTICE[2] );
    glVertex3d( LATTICE[3] , LATTICE[4] , LATTICE[5] );
    glVertex3d( LATTICE[3] + LATTICE[6] - LATTICE[0] ,
                LATTICE[4] + LATTICE[7] - LATTICE[1] ,
                LATTICE[5] + LATTICE[8] - LATTICE[2] );
    glVertex3d( LATTICE[6] , LATTICE[7] , LATTICE[8] );
    glVertex3d( LATTICE[0] , LATTICE[1] , LATTICE[2] );
    glVertex3d( LATTICE[9] , LATTICE[10], LATTICE[11] );
    glVertex3d( LATTICE[6] + LATTICE[9] - LATTICE[0] ,
                LATTICE[7] + LATTICE[10]- LATTICE[1] ,
                LATTICE[8] + LATTICE[11]- LATTICE[2] );
    glVertex3d( LATTICE[6] , LATTICE[7] , LATTICE[8] );
    glVertex3d( LATTICE[6] + LATTICE[9] - LATTICE[0] ,
                LATTICE[7] + LATTICE[10]- LATTICE[1] ,
                LATTICE[8] + LATTICE[11]- LATTICE[2] );
    glVertex3d( LATTICE[3] + LATTICE[6] + LATTICE[9]-2.0*LATTICE[0],
                LATTICE[4] + LATTICE[7] + LATTICE[10]-2.0*LATTICE[1],
                LATTICE[5] + LATTICE[8] + LATTICE[11]-2.0*LATTICE[2]);
    glVertex3d( LATTICE[3] + LATTICE[6] - LATTICE[0] ,
                LATTICE[4] + LATTICE[7] - LATTICE[1] ,
                LATTICE[5] + LATTICE[8] - LATTICE[2] );
    glVertex3d( LATTICE[3] + LATTICE[6] + LATTICE[9]-2.0*LATTICE[0] ,
                LATTICE[4] + LATTICE[7] + LATTICE[10]-2.0*LATTICE[1],
                LATTICE[5] + LATTICE[8] + LATTICE[11]-2.0*LATTICE[2]);
    glVertex3d( LATTICE[3] + LATTICE[9] - LATTICE[0] ,
                LATTICE[4] + LATTICE[10]- LATTICE[1] ,
                LATTICE[5] + LATTICE[11]- LATTICE[2] );
    glVertex3d( LATTICE[3] , LATTICE[4] , LATTICE[5] );
    glVertex3d( LATTICE[3] + LATTICE[9] - LATTICE[0] ,
                LATTICE[4] + LATTICE[10]- LATTICE[1] ,
                LATTICE[5] + LATTICE[11]- LATTICE[2] );
    glVertex3d( LATTICE[9] , LATTICE[10], LATTICE[11] );
    glEnable( GL_LIGHTING );
  }
  else if (surface[n]->type==CURVE) // AQUI
  {

     for ( i = mini+1 ; i < maxi ; i++ )
     {
        grid->GetCubePoints(i,minj+1,mink+1,grdvec);

        surface[n]->GetCubeValues(i,minj+1,mink+1,grdval);

        Vector v1 = grdvec[2];
        v1 -= grdvec[3];
        Vector v2 = grdvec[1];
        v2 -= grdvec[2];
        v2 <= v1;

        v1 = v2;
        v1 ^= grdval[3]*5;

        v1 += grdvec[3];
        glVertex3d(v1[0],v1[1],v1[2]);

        v2 ^= grdval[2]*5;

        v2 += grdvec[2];
        glVertex3d(v2[0],v2[1],v2[2]);
	v1[2] = BohrToAngstrom(v1[2]);

	cout << v1[2] << "	" << grdval[3] << endl;
     }
  }
  else
  {
    adjust_max = -1.0/GN_TOL_NUM;
    adjust_min = -adjust_max;
    for ( k = mink+1 ; k < maxk ; k++ )
    {
      for ( j = minj+1 ; j < maxj ; j++ )
      {
        for ( i = mini+1 ; i < maxi ; i++ )
        {
          // Get grid points for a cube, square or line
          grid->GetCubePoints(i,j,k,grdvec);
          surface[n]->GetCubeValues(i,j,k,grdval);
          // Search for triangles
          if (surface[n]->type==COLOURED_PLANE||
              surface[n]->type==CONTOUR)
          {
            if (( plane_axis == 0 && plane_point == i )||
                ( plane_axis == 1 && plane_point == j )||
                ( plane_axis == 2 && plane_point == k )) 
            {
              if (surface[n]->type==COLOURED_PLANE)
              {
                int ii;
                int ptr[3][4] = {{0,3,7,4},{0,1,5,4},{0,1,2,3}};
                // glNormal3d(dstep[2][0],dstep[2][1],dstep[2][2]);
                ii = ptr[plane_axis][0];
                RealToRGB(grdval[ii],rgb);
                glColor3dv(rgb);
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 0

                ii = ptr[plane_axis][1];
                RealToRGB(grdval[ii],rgb);
                glColor3dv(rgb);
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 1

                ii = ptr[plane_axis][3];
                RealToRGB(grdval[ii],rgb);
                glColor3dv(rgb);
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 3
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 3

                ii = ptr[plane_axis][1];
                RealToRGB(grdval[ii],rgb);
                glColor3dv(rgb);
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 1

                ii = ptr[plane_axis][2];
                RealToRGB(grdval[ii],rgb);
                glColor3dv(rgb);
                glVertex3d(grdvec[ii][0],grdvec[ii][1],grdvec[ii][2]);  // 2
              }
              else if (surface[n]->type==CONTOUR)
              {
                int nl;
                double lines[2][2][3];
                LBox(surface[n]->iso,lines,&nl);
                for (int il = 0; il < nl ; il++)
                {
                  glVertex3dv(&lines[il][0][0]);
                  glVertex3dv(&lines[il][1][0]);
                }
              }
            }
            continue;
          }
          else TBox(surface[n]->iso,t,&nt);
          if (nt>0) 
          {
            int nv = 0;
            double x[15],y[15],z[15],v[15];

            if (surface[n]->type==COLOURED_ISO||
                surface[n]->type==COLOURED_PLANE)
            {
              nv = 3*nt;
              int iv=0;
              for ( it = 0 ; it < nt ; it++ )
              {
                for ( ip = 0 ; ip < 3 ; ip++ )
                {
                  x[iv] = t[ip][it][0];
                  y[iv] = t[ip][it][1];
                  z[iv] = t[ip][it][2];
                  iv++;
                }
              }
              surface[n]->GetValues(x,y,z,v,nv);
            }

            // Got some triangles!
            GetNormals(n,t,tn,nt);
            jgp = 0;
            for ( it = 0 ; it < nt ; it++ )
            {
              for ( ip = 0 ; ip < 3 ; ip++ )
              {
                if (surface[n]->type==COLOURED_ISO)
                {
                  RealToRGB(v[jgp+ip],rgb);
                  adjust_max = GN_MAX(v[jgp+ip],adjust_max);
                  adjust_min = GN_MIN(v[jgp+ip],adjust_min);
                  glColor3dv(rgb);
                }
                glNormal3dv(&tn[ip][it][0]);
                glVertex3dv(&t[ip][it][0]);
              }
              jgp = jgp + 3;
            }
          }
        }
      }
    }
  }

}

// Purpose: Get isosurface triangles from a cubic box.
//
// (c) This an adaption of Paul Bourke's code.
void Plotter::TBox(double iso, double t[3][5][3], int* nt)
{
//
//    Order of the cube vertices:
//
//         4----------5    
//        /|         /|   
//       / |        / |
//      7----------6  |  
//      |  |       |  | 
//      |  |       |  |   
//      |  0-------|--1  
//      | /        | /  
//      |/         |/  
//      3----------2
//

  int i,j,k,indx; 
  double a,b;
  double vlist[12][3];

  int edgetable[256]={
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

  int tritable[256][16] =
  {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  int ptr[2][12] = {{0,1,2,3,4,5,6,7,0,1,2,3},
                    {1,2,3,0,5,6,7,4,4,5,6,7}};

  *nt = 0;

  indx = 0;
  j = 1;
  for ( i = 0 ; i < 8 ; i++ )
  {
    if (grdval[i]<iso) indx = indx + j;
    j = j + j;
  }

  if (indx==0||indx==255) return;

  j = 1;
  for ( i = 0 ; i < 12 ; i++ )
  {
    if (edgetable[indx] & j)
    {
      a = grdval[ptr[0][i]];
      b = grdval[ptr[1][i]];
      // Get our approximate intersection ISO values
      if (GN_ABS(iso-a)<GN_TOL_NUM) 
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else if(GN_ABS(iso-b)<GN_TOL_NUM) 
      {
        vlist[i][0] = grdvec[ptr[1][i]][0];
        vlist[i][1] = grdvec[ptr[1][i]][1];
        vlist[i][2] = grdvec[ptr[1][i]][2];
      }
      else if (GN_ABS(b-a)<GN_TOL_NUM)
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else
      {
        vlist[i][0] = grdvec[ptr[0][i]][0] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][0]-grdvec[ptr[0][i]][0]);
        vlist[i][1] = grdvec[ptr[0][i]][1] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][1]-grdvec[ptr[0][i]][1]);
        vlist[i][2] = grdvec[ptr[0][i]][2] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][2]-grdvec[ptr[0][i]][2]);
      }
    }
    j = j + j;
  }

  i = 0; 
  while ( tritable[indx][i] != -1)
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      for ( k = 0 ; k < 3 ; k++ )
      {
        t[j][*nt][k] = vlist[tritable[indx][i+j]][k];
      }
    }
    (*nt)++;
    i += 3; 
  }
}

void Plotter::GetNormals(int n,double rt[3][5][3],double rtn[3][5][3],int rnt)
{
  int i,rit,riv;
  double x[15];
  double y[15];
  double z[15];
  double dx[15];
  double dy[15];
  double dz[15];
  int nr;

  nr = 0;
  for (rit=0;rit<rnt;rit++)
  {
    for (riv=0;riv<3;riv++)
    {
      x[nr] = rt[riv][rit][0];
      y[nr] = rt[riv][rit][1];
      z[nr] = rt[riv][rit][2];
      nr++;
    }
  }
  surface[n]->GetDerivatives(x,y,z,dx,dy,dz,nr);
  i = 0;
  for (rit=0;rit<rnt;rit++)
  {
    for (riv=0;riv<3;riv++)
    {
      rtn[riv][rit][0] = dx[i]*surface[n]->iso;
      rtn[riv][rit][1] = dy[i]*surface[n]->iso;
      rtn[riv][rit][2] = dz[i]*surface[n]->iso;
      i++;
    }
  }
}

// Get contour lines in a square box.
void Plotter::LBox(double iso,double l[2][2][3],int* nl)
{
//
// 
//
//            0 ----------- 1
//           /             /
//          /             /
//         /             /
//        3 ----------- 2
//
  int i,il,ip,j,indx;
  int couple[2][2][2];
  double EPS=0.00001;

  indx = 0;
  j = 1;
  for (i=0;i<4;i++)
  {
    if (grdval[i]<iso) indx += j;
    j *= 2;
  }
  if (indx>7) indx = 15 - indx;

  *nl = 1;
  if (indx==1) 
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 0;
    couple[0][1][1] = 3;
  }
  else if (indx==2)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
  }
  else if (indx==3)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 3;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
  }
  else if (indx==4) 
  {
    couple[0][0][0] = 1;
    couple[0][0][1] = 2;
    couple[0][1][0] = 2;
    couple[0][1][1] = 3;
  }
  else if (indx==5)
  {
    *nl = 2;
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
    couple[1][0][0] = 2;
    couple[1][0][1] = 3;
    couple[1][1][0] = 3;
    couple[1][1][1] = 0;
  }
  else if (indx==6)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 2;
    couple[0][1][1] = 3;
  }
  else if (indx==7)
  {
    couple[0][0][0] = 2;
    couple[0][0][1] = 3;
    couple[0][1][0] = 3;
    couple[0][1][1] = 0;
  }
  else if (indx==8)
  {
    *nl = 2;
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 0;
    couple[0][1][1] = 3;
    couple[1][0][0] = 1;
    couple[1][0][1] = 2;
    couple[1][1][0] = 2;
    couple[1][1][1] = 3;
  }
  else *nl = 0;
        
  double factor;
  for ( il = 0 ; il < *nl ; il++ )
  {
    for ( ip = 0 ; ip < 2 ; ip++ )
    {
      i = couple[il][ip][0];
      j = couple[il][ip][1];
      if (GN_ABS(iso-grdval[i])<EPS)
      {
        l[il][ip][0] = grdvec[i][0];
        l[il][ip][1] = grdvec[i][1];
        l[il][ip][2] = grdvec[i][2];
      }
      else if (GN_ABS(iso-grdval[j])<EPS)
      {
        l[il][ip][0] = grdvec[j][0];
        l[il][ip][1] = grdvec[j][1];
        l[il][ip][2] = grdvec[j][2];
      }
      else if (GN_ABS(grdval[i]-grdval[j])<EPS)
      {
        l[il][ip][0] = grdvec[i][0];
        l[il][ip][1] = grdvec[i][1];
        l[il][ip][2] = grdvec[i][2];
      }
      else 
      {
        factor = (iso-grdval[i])/(grdval[j]-grdval[i]);
        l[il][ip][0] = grdvec[i][0] + factor*(grdvec[j][0]-grdvec[i][0]);
        l[il][ip][1] = grdvec[i][1] + factor*(grdvec[j][1]-grdvec[i][1]);
        l[il][ip][2] = grdvec[i][2] + factor*(grdvec[j][2]-grdvec[i][2]);
      }
    }
  }
}

// Translate real value (0,1) to RGB color.
void Plotter::RealToRGB(double r,double *rgb)
{
  int n,irgb[3];
  double eps = 0.00000001;

  n = int(255.0*GN_MIN(1.0,(r-fmin)/(fmax-fmin+eps)));
  
  if (n<128) irgb[0] = 0;
  else if (n<192) irgb[0] = (n-128)*4 + 1;
  else irgb[0] = 255;
 
  if (n<1) irgb[1] = 0;
  else if (n<64) irgb[1] = 4*n-1;
  else if (n<192) irgb[1] = 255;
  else irgb[1] = 255 - 4*(n-192) - GN_MOD(n,2);

  if (n<64) irgb[2] = 255;
  else if (n<128) irgb[2] = 255-4*(n-64)-GN_MOD(n,2);
  else irgb[2] = 0;

  rgb[0] = irgb[0]/255.0;
  rgb[1] = irgb[1]/255.0;
  rgb[2] = irgb[2]/255.0;
}



void Plotter::Clear()
{
  surface.clear();
}

void Plotter::ChangeStyle( const QString &newstyle )
{
  for ( int i = 0 ; i < NPLOT_STYLES ; i++ )
  {
    if ( newstyle == stylenames[i] )
    {
      style = i;
      break;
    };
  };
}

void Plotter::ChangeType( const QString &newtype )
{
  for ( int i = 0 ; i < NPLOT_TYPES ; i++ )
  {
    if ( newtype == typenames[i] )
    {
      type = i;
      break;
    };
  };
}

void Plotter::ChangePlaneAxis( int npa )
{
  plane_axis = npa-1;
}

void Plotter::ChangePlanePoint( int npp )
{
  plane_point = npp-1;
}

void Plotter::ChangeContourMin( double cmin ) 
{ fmin = cmin; }

void Plotter::ChangeContourMax( double cmax ) 
{ fmax = cmax; }

void Plotter::ChangeContourStep( double cstep ) 
{ fstep = cstep; }


void Plotter::Adjust(double *iso)
{
  int n = surface.size() - 1;
  double smin,smax;
  surface[n]->GetExtrema(&smin,&smax);
  fmin = smin;
  fmax = smax;
  fstep = (fmax - fmin)/30.0;
  *iso = (fmax + fmin)/2.0;
}


