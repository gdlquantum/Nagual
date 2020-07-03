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
#ifndef GN_DWPANEL_H
#define GN_DWPANEL_H

#include <vector>

#include <qgl.h>
#include <qmainwindow.h>

using namespace std;

class QSpinBox;
class QComboBox;
class QToolButton;
class QDoubleSpinBox;
class QResizeEvent;
class QCheckBox;
class QTabWidget;

class GNagual; 
class System;
class Viewer;
class Model;
class DwGeo;
class Plotter;
class DwPlot;

class DwPanel : public QWidget
{
  Q_OBJECT

  public:

    DwPanel( GNagual*, System* );

    void Report( const char*, ... );

    GNagual *gn;
    System *sys;
    Viewer *viewer;
    DwGeo *geoeditor;
    Model* model;
    QMenuBar* menu;
    DwPlot *winplot;

  public slots:
 
    void ChangeModel( const QString & );
    void About( void );
    void MakeStyle( const QString & );
    void ReadFile( void );
    void DoReadFile(char*,char*);
    void WriteFile( void );
    void DoWriteFile(char*,char*);

  protected:

    void CreateActions();
    void CreateMenus();

    QStatusBar *statusbar;

    QMenuBar *mainMenu;

    QMenu *fileMenu;
    QMenu *colorMenu;
    QMenu *readMenu;
    QMenu *writeMenu;
    QMenu *helpMenu;

    vector<QAction*> RFMTA;
    QAction *quitAct;
    QAction *aboutAct;
    QAction *imgAct;
    QAction *coloratomAct;
    QAction *colorbgAct;
    QAction *colortagAct;
    QAction *colorarrowAct;
    QAction *colorsurfAct;
    QAction *exitAction;

    QToolButton *tstButton;
    QToolButton *inButton;
    QToolButton *outButton;
    QToolButton *giraButton;
    QToolButton *numAtomButton;

    QComboBox *modelComboBox;
    QComboBox *styleComboBox;
    QComboBox *arrowComboBox;
   
    QTimer* timer;
    QSpinBox *lineS;

    QTabWidget *tabWidget;

  protected slots:

    void ChangeAtomColor(void);
    void ChangeSurfaceColor(void);
    void WriteImage(void);
    void SetSpin(void);
    void AdvanceSpin(void);

};

#endif // GN_DWPANEL_H

