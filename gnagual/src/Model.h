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
#ifndef GN_MODEL_H
#define GN_MODEL_H

#include <GL/glu.h>

#include <QtCore>

class DwPanel;
class Atom;
class Molecule;

class Model
{
  public:

    Model(DwPanel*);

    DwPanel* panel;
    Molecule *mol;
    double srs;
    double crs;
    bool drawmol;
    GLUquadricObj* q;

    void DrawSphere( double* , double, int );
    void DrawCylinder( double* , double* , double, int );
    void DrawCone( double* , double* , double, int );
    void DrawMolecule( int, bool );
    void ChangeKind( const QString & );

  public slots:

    void ChangeArrows( const QString & );
    void ChangeArrowColor( void );

  protected:

    QString arrow_type;
    void DrawAtom( Atom* , double, int );
    void DrawBond( Atom* , Atom* , double, int );
};

#endif // GN_MODEL_H
