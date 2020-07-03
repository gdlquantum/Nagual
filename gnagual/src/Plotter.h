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
#ifndef GN_PLOTTER_H
#define GN_PLOTTER_H

#include <QtGui>

#include <vector>

#include <Vector.h>

#define ISO             0
#define CONTOUR         1
#define COLOURED_ISO    2
#define COLOURED_PLANE  3
#define CURVE           4
#define NPLOT_TYPES     5

#define SOLID_STYLE         0
#define POINTS_STYLE        1
#define TRANSLUCENT_STYLE   2
#define MESH_STYLE          3
#define NPLOT_STYLES        4

extern char* stylenames[NPLOT_STYLES];
extern char* typenames[NPLOT_TYPES];

using namespace std;

class GNagual;
class DwPlot;
class Grid;
class Surface;

class Plotter : public QWidget
{
  Q_OBJECT

  public:
    Plotter(GNagual*);

    GNagual *gn;
    Grid *grid;
    vector<Surface*> surface;
    double fmin;
    double fmax;
    double fstep;
    double adjust_max;
    double adjust_min;
    bool use_lists;

    void Adjust(double*);
    void Clear(void);
    void GetColor( int , double* );
    void DoSurface( int );
    void SetColor( int , double* );
    int NSurf( void ); 

  public slots:

    void ChangeStyle( const QString & );
    void ChangeType( const QString & );
    void ChangePlaneAxis( int );
    void ChangePlanePoint( int );
    void ChangeContourMin( double );
    void ChangeContourMax( double );
    void ChangeContourStep( double );
    void RequestPlot(const QString&,int,int,int,double,bool);  

  protected:


    void Build(int);
    void WriteCubeFile(int);
    void PlotISO(int);
    void PlotCURVE(int);
    void PlotCOLOURED_ISO(int);
    void PlotCOLOURED_PLANE(int);
    void PlotCONTOUR(int);
    void TBox(double iso,double t[3][5][3],int* nt);
    void LBox(double iso,double l[2][2][3],int* nl);
    void GetNormals(int,double t[3][5][3], double tn[3][5][3],int);
    void RealToRGB(double,double*);

    int line_width;
    int point_size;
    int shiny;
    int style;
    int type;

    int plane_axis;
    int plane_point;

    double grdval[8];
    Vector grdvec[8];
};

#endif  // GN_PLOTTER_H

