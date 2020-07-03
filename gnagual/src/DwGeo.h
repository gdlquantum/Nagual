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
#ifndef GN_DWGEO_H
#define GN_DWGEO_H

#include <vector>

#include <QtGui>
#include <QStringList>

#include <Parameter.h>
//Qt5 port
#include <QTableWidgetItem>

#define Cartesians 0
#define ZMatrix    1

using namespace std;

class GNagual;
class DwPanel;
class Molecule;

class DwGeo : public QWidget
{
  Q_OBJECT

  public:
    DwGeo( DwPanel* );

    DwPanel *panel;
    GNagual *gn;
    Molecule *mol;
    int default_element;
    bool use_default;
    char gname[80];

    void Enabled( bool );
    void AddAtom( int, int, int );
    void AddGroup(int,int,int);
    void Read(char*,char*);

  public slots:

    void AddAtomStart( void );
    void AddGroupStart(const QString &);
    void DeleteAtom( void );
    void Edited(QTableWidgetItem*);
    void SetupMeasurement(const QString &);
    void BuildGroup( const QString & );    
    void ChangeCoordinates( void );
    void ChangeUnits( void );
    void Update( void );
    void FitDistance( void );
    void UseDefault( void );
    int SGDrv(int*);

  protected:

    bool bohr;
    char group[GN_MAX_STR_SIZE];
    int CCSET;
    QTableWidget *tw;
    QStringList *Headers;
    QLineEdit *smiles;
    QTimer* timer;
};

#endif // GN_DWGEO_H

