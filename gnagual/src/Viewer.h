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
#ifndef GN_VIEWER_H
#define GN_VIEWER_H

#include <vector>
#include <QWheelEvent>
#include <QtOpenGL/QtOpenGL>

#include <Vector.h>

using namespace std;

class wheelEvent;

class GNagual;
class DwPanel;
class Plotter;
class DwPlot;
class DwImage;

class Viewer : public QGLWidget
{
  Q_OBJECT

  public:
    Viewer( GNagual*, DwPanel* );

    GNagual *gn;
    DwPanel *panel;
    Plotter* plotter;
    DwImage *render;

    double trans[3];
    double scale;
    GLdouble rot[16];

    void Project( void );
    void PickSurface( void );

    void SaveSelection( int );
    void DrawSurfaces( void );
    void DrawStrings( void );
    void Rotate( double , double , double , double );
    void ResetRotation(void);
    void Write(char*,Vector,size_t);

    vector<vector<size_t> > monitored_bonds;
    vector<vector<size_t> > monitored_angles;
    vector<vector<size_t> > monitored_dihedrals;

  public slots:
 
    void Redraw( void );
    void ZoomIn( void );
    void ZoomOut( void );
    void ChangeBackgroundColor( void );
    void ChangeTags( const QString & );
    void ChangeTagColor( void );
    void wheelEvent(QWheelEvent *event);

  protected:

    QString tag_type;
    QPrinter *printer;
    int mouse[2];
    float view_width;
    float view_height;
    float view_far;
    float view_near;

    void initializeGL( void );
    void resizeGL( int , int );
    void paintGL( void );

    int SelectedAtom( size_t );

    double tag_color[4];
    double arrow_color[4];
    double bg[4];
    vector<size_t> selected_atoms;

  protected slots:

    void mousePressEvent( QMouseEvent* );
    void mouseMoveEvent( QMouseEvent* );
    void mouseReleaseEvent( QMouseEvent* );
    void keyPressEvent( QKeyEvent* );

    void PickAtom( void );

};

#endif  // GN_VIEWER_H

