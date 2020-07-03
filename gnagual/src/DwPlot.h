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
#ifndef GN_DWPLOT_H
#define GN_DWPLOT_H

#include <vector>

#include <QWidget>

#include <Vector.h>

using namespace std;

class QListWidget;
class QTabWidget;
class QVBoxLayout;
class QDoubleSpinBox;
class QCheckBox;
class QRadioButton;

class GNagual;
class Viewer;
class Vector;
class Plotter;

class DwPlot : public QWidget
{
  Q_OBJECT

  public:
    DwPlot( GNagual* , Viewer* );

    GNagual *gn;
    Viewer *viewer;
    Plotter *plotter;

    vector<double> temp;

  public slots:

    void Adjust(void);
    void ChangeIso( double );
    void ChangeLineWidth( int );
    void ChangePointSize( int );
    void ChangeBox( double );
    void ChangeMesh( double );
    void ChangeFlags( void );  
    void ChangeSetNumber( int );
    void ChangeOrbitalNumber( int );
    void ChangePerturbationNumber( int );
    void SetupPlot(void);
    void BuildPlot(void);

  protected:

    void NewList( char** , int , char* ); 
    void SetupBox( QVBoxLayout*, double );

    int hold;
    int hide_surfaces;
    int drawref;
    int approximate_normals;
    int average_normals;

    int shiny;
    int line_width;
    int set_number;
    int orbital_number;
    int perturbation_number;
    int point_size;
    int lastplaneaxis;
    double isovalue;

    QDoubleSpinBox *cminS;
    QDoubleSpinBox *cmaxS;
    QDoubleSpinBox *cstepS;
    QDoubleSpinBox *isoS;
    QDoubleSpinBox *meshS;

    QDoubleSpinBox* axis[4][3];
    QCheckBox* flagbox[4][3];
    QTabWidget *tabWidget;
    vector<QWidget*> tabs;
    vector<QListWidget*> listados;
    vector<double> grdval;
    vector<double> grdvaldegen;
    vector<Vector> grdvec;
    vector<Vector> grdvecdegen;

};

#endif  // GN_DWPLOT_H

