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
#include <DwPlot.h>
#include <Viewer.h>
#include <DwPanel.h>
#include <GNagual.h>
#include <Vector.h>
#include <Plotter.h>
#include <Grid.h>
#include <Surface.h>
#include <System.h>

using namespace std;

static char *flagbox_names[4][3] = {
{"Shiny","Reference","Nothing"},
{"Hide all","Hold Previous","Nothing"},
{"Approximate","Average","Nothing"},
{"Nothing","Nothing","Nothing"}};

DwPlot::DwPlot( GNagual *ix , Viewer* view)
       : QWidget( 0 )
{
  gn = ix;
  viewer = view;
  plotter = viewer->plotter;
  tabWidget = new QTabWidget;

  lastplaneaxis = 2;

  QGroupBox *settBox = new QGroupBox(tr("Settings"));
  QHBoxLayout *settLayout = new QHBoxLayout;
  settLayout->addStretch(1);

  QComboBox* typeC = new QComboBox( settBox );
  for ( int i = 0 ; i < NPLOT_TYPES ; i++ )
  {
    typeC->addItem( typenames[i] );
  }; 
  connect( typeC , SIGNAL( activated( const QString & ) ) , 
           plotter , SLOT( ChangeType( const QString & ) ) );
  settLayout->addWidget( typeC );

  QComboBox* styleC = new QComboBox( settBox );
  for ( int i = 0 ; i < NPLOT_STYLES ; i++ )
  {
    styleC->addItem( stylenames[i] );
  }; 
  connect( styleC , SIGNAL( activated( const QString & ) ) , 
           plotter , SLOT( ChangeStyle( const QString & ) ) );
  settLayout->addWidget( styleC );

  QPushButton *freshB = new QPushButton(tr("Adjust"));
  connect( freshB, SIGNAL(clicked()), this, SLOT(Adjust()));
  settLayout->addWidget( freshB );

  QPushButton *plotB = new QPushButton(tr("Build"));
  connect( plotB, SIGNAL(clicked()), this, SLOT(BuildPlot()));
  settLayout->addWidget( plotB );

  QPushButton *setB = new QPushButton(tr("Compute"));
  connect( setB, SIGNAL(clicked()), this, SLOT(SetupPlot()));
  settLayout->addWidget( setB );

  settBox->setLayout( settLayout );

  QGroupBox *valBox = new QGroupBox(tr("Values: ISO/MIN/MAX/STEP"));
  QHBoxLayout *valLayout = new QHBoxLayout;
  valLayout->addStretch(1);

  double valrange = 10000.0;
  isovalue = 0.02;
  isoS = new QDoubleSpinBox( valBox );
  isoS->setRange( -valrange , valrange );
  isoS->setSingleStep( 0.00001 );
  isoS->setDecimals( 10 );
  isoS->setValue( isovalue );
  connect( isoS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeIso( double ) ) );
  valLayout->addWidget( isoS );

  cminS = new QDoubleSpinBox( valBox );
  cminS->setRange( -valrange , valrange );
  cminS->setSingleStep( 0.01 );
  cminS->setDecimals( 10 );
  cminS->setValue( -1.0 );
  valLayout->addWidget( cminS );
  connect( cminS , SIGNAL(valueChanged( double )) ,
           plotter , SLOT( ChangeContourMin( double ) ) );

  cmaxS = new QDoubleSpinBox( valBox );
  cmaxS->setRange( -valrange , valrange );
  cmaxS->setSingleStep( 0.01 );
  cmaxS->setValue( 1.0 );
  cmaxS->setDecimals( 10 );
  valLayout->addWidget( cmaxS );
  connect( cmaxS , SIGNAL(valueChanged( double )) ,
           plotter , SLOT( ChangeContourMax( double ) ) );

  cstepS = new QDoubleSpinBox( valBox );
  cstepS->setRange( -valrange , valrange );
  cstepS->setSingleStep( 0.001 );
  cstepS->setDecimals( 10 );
  cstepS->setValue( 0.1 );
  valLayout->addWidget( cstepS );
  connect( cstepS , SIGNAL(valueChanged( double )) ,
           plotter , SLOT( ChangeContourStep( double ) ) );

  valBox->setLayout(valLayout);

  QGroupBox *orbBox = new QGroupBox(tr("Set/Number/Perturbation"));
  QHBoxLayout *orbLayout = new QHBoxLayout;
  orbLayout->addStretch(1);

  set_number = 0;
  QSpinBox *setS = new QSpinBox( orbBox );
  setS->setRange( 0 , GN_MAX_NSPIN - 1 );
  setS->setSingleStep( 1 );
  setS->setValue( set_number );
  orbLayout->addWidget( setS );
  connect( setS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeSetNumber( int ) ) );

  orbital_number = 0;
  QSpinBox *onS = new QSpinBox( orbBox );
  onS->setRange( 0 , 1000 );
  onS->setSingleStep( 1 );
  onS->setValue( orbital_number );
  orbLayout->addWidget( onS );
  connect( onS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeOrbitalNumber( int ) ) );

  perturbation_number = 0;
  QSpinBox *pnS = new QSpinBox( orbBox );
  pnS->setRange( 0 , 1000 );
  pnS->setSingleStep( 1 );
  pnS->setValue( perturbation_number );
  orbLayout->addWidget( pnS );
  connect( pnS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangePerturbationNumber( int ) ) );

  orbBox->setLayout( orbLayout );

  QGroupBox *drawBox = new QGroupBox(tr("Drawing: Grid/Line/Point"));
  QHBoxLayout *drawLayout = new QHBoxLayout;
  drawLayout->addStretch(1);

  meshS = new QDoubleSpinBox( drawBox );
  meshS->setRange( 0.02*GN_GS , 2*GN_GS );
  meshS->setSingleStep( 0.02*GN_GS );
  meshS->setDecimals( 5 );
  meshS->setValue( 0.1*GN_GS );
  drawLayout->addWidget( meshS );
  connect( meshS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeMesh( double ) ) );

  line_width = 2;
  QSpinBox *lineS = new QSpinBox( drawBox );
  lineS->setRange( 1 , 10 );
  lineS->setSingleStep( 1 );
  lineS->setValue( line_width );
  drawLayout->addWidget( lineS );
  connect( lineS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeLineWidth( int ) ) );

  QSpinBox *psS = new QSpinBox( drawBox );
  psS->setRange( 1 , 10 );
  psS->setSingleStep( 1 );
  psS->setValue( 3 );
  drawLayout->addWidget( psS );
  connect( lineS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangePointSize( int ) ) );

  drawBox->setLayout( drawLayout );

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget( tabWidget );
  mainLayout->addWidget( orbBox );
  mainLayout->addWidget( drawBox );
  mainLayout->addWidget( valBox );
  mainLayout->addWidget( settBox );

  SetupBox( mainLayout , GN_GS );
  setLayout( mainLayout );
  setWindowTitle(tr( "Plot Controls" ));

  tabs.clear();
  listados.clear();
  tabWidget->clear();

  int imo;
  int NMAX = 50;
  char *names[NMAX];

  for ( imo = 0 ; imo  < NMAX ; imo++ )
    names[imo] = (char*) malloc(256*sizeof(char));

  int ns = 0; 
  sprintf( names[ns],"BOX"); ns++;
  sprintf( names[ns],"GTO"); ns++;
  sprintf( names[ns],"ORBITAL"); ns++;
  sprintf( names[ns],"DENSITY"); ns++;
  sprintf( names[ns],"CHARGE DENSITY"); ns++;
  sprintf( names[ns],"SPECIES DENSITY"); ns++;
  sprintf( names[ns],"SPIN DENSITY"); ns++;
  sprintf( names[ns],"SHANON ENTROPY"); ns++;

  NewList( names , ns , "" );

  for ( imo = 0 ; imo  < NMAX ; imo++ )
    free(names[imo]);
}

void DwPlot::NewList( char** entries , int nentries , char* name ) 
{
  tabs.push_back( new QWidget( 0 ) );
  QListWidget *slListBox = new QListWidget;
  QStringList sllist;

  for ( int i = 0; i < nentries; i++ )
  {
    sllist.append(tr(entries[i]));
  };
  slListBox->insertItems( 0, sllist);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget( slListBox );
  tabs[tabs.size()-1]->setLayout( layout );
  listados.push_back( slListBox );
  tabWidget->addTab( tabs[tabs.size()-1] , QString(name) );
}

void DwPlot::ChangeMesh( double newval ) 
{ ChangeBox( newval ); }

void DwPlot::ChangeIso( double newiso ) 
{ isovalue = newiso; }

void DwPlot::ChangeLineWidth( int lw ) 
{ line_width = lw; }

void DwPlot::ChangePointSize( int ps ) 
{ point_size = ps; }





void DwPlot::SetupBox( QVBoxLayout *mainLayout, double sf )
{
  int i,j,k;

  QTabWidget *boxtab = new QTabWidget;
  QWidget *cube = new QWidget( 0 );
  QVBoxLayout *cubeLayout = new QVBoxLayout;

  double axisval[] = {-sf,-sf,-sf,
                       sf,-sf,-sf,
                      -sf, sf,-sf,
                      -sf,-sf, sf};

  char *cube_names[4] = {"Origin", "Axis 1", "Axis 2", "Axis 3"};
  QGroupBox* cube_g[4];
  QHBoxLayout* cube_l[4]; 

  k = 0;
  for ( i = 0 ; i < 4; i++ )
  {
    cube_g[i] = new QGroupBox(tr(cube_names[i]));
    cube_l[i] = new QHBoxLayout;
    cube_l[i]->addStretch(1);
    for ( j = 0 ; j < 3 ; j++ )
    {
      axis[i][j] = new QDoubleSpinBox( cube_g[i] );
      axis[i][j]->setRange( -100000.0 , 100000.0 );
      axis[i][j]->setSingleStep( 1.0 );
      axis[i][j]->setValue( axisval[k] );
      axis[i][j]->setDecimals( 5 );
      connect( axis[i][j] , SIGNAL(valueChanged( double )) ,
               this , SLOT( ChangeBox( double ) ) );
      cube_l[i]->addWidget( axis[i][j] );
      k++;
    };
    cube_g[i]->setLayout( cube_l[i] );
    cubeLayout->addWidget( cube_g[i] );
  };
  cube->setLayout( cubeLayout );
  boxtab->addTab( cube , tr("Box") );

  QWidget *flags = new QWidget( 0 );
  QVBoxLayout *flagsLayout = new QVBoxLayout;

  char *flags_names[4] = {"General","Visible","Normals","Nothing"};
  QGroupBox* flags_g[4];
  QHBoxLayout* flags_l[4]; 

  for ( i = 0 ; i < 4; i++ )
  {
    flags_g[i] = new QGroupBox(tr( flags_names[i]));
    flags_l[i] = new QHBoxLayout;
    flags_l[i]->addStretch(1);
    for ( j = 0 ; j < 3; j++ )
    {
      flagbox[i][j] = new QCheckBox( flags_g[i] );
      flagbox[i][j]->setText( flagbox_names[i][j] ); 
      connect( flagbox[i][j] , SIGNAL( clicked() ) , 
               this , SLOT( ChangeFlags() ) );
      flags_l[i]->addWidget( flagbox[i][j] );
    };
    flags_g[i]->setLayout( flags_l[i] );
    flagsLayout->addWidget( flags_g[i] );
  };
  flags->setLayout( flagsLayout );
  boxtab->addTab( flags , tr("Flags") );

  QWidget *plane = new QWidget( 0 );
  QVBoxLayout *planeLayout = new QVBoxLayout;

  QGroupBox *axisBox = new QGroupBox(tr("Axis"));
  QHBoxLayout *axisLayout = new QHBoxLayout;
  axisLayout->addStretch(1);

   QSpinBox *axisS = new QSpinBox( axisBox );
   axisS->setRange( 1 , 3 );
   axisS->setSingleStep( 1 );
   axisS->setValue( 3 );
   connect( axisS , SIGNAL(valueChanged( int )) ,
            plotter , SLOT( ChangePlaneAxis( int ) ) );


  axisLayout->addWidget( axisS );
  axisBox->setLayout( axisLayout );

  QGroupBox *pointBox = new QGroupBox(tr("Point"));
  QHBoxLayout *pointLayout = new QHBoxLayout;
  pointLayout->addStretch(1);

  QSpinBox *pointS = new QSpinBox( pointBox );
  pointS->setRange( 1 , 100000 );
  pointS->setSingleStep( 1 );
  pointS->setValue( 1 );
  connect( pointS , SIGNAL(valueChanged( int )) ,
           plotter , SLOT( ChangePlanePoint( int ) ) );

  pointLayout->addWidget( pointS );
  pointBox->setLayout( pointLayout );

  planeLayout->addWidget( axisBox );
  planeLayout->addWidget( pointBox );
  plane->setLayout( planeLayout );
  boxtab->addTab( plane , tr("Plane") );

  mainLayout->addWidget( boxtab );
  hold = 0;
  drawref = 0;
  approximate_normals = 0;
  average_normals = 0;
  hide_surfaces = 0;
  shiny = 0;

  plotter->grid->SetLattice(axisval,meshS->value());
}

void DwPlot::ChangeBox( double val )
{
  int i,j;

  double lat[12];
  for ( i = 0 ; i < 4; i++ )
  {
    for ( j = 0 ; j < 3; j++ )
    {
      lat[3*i+j] = axis[i][j]->value();
    };
  };
  plotter->grid->SetLattice(lat,meshS->value());
}

void DwPlot::ChangeFlags()
{
  int i,j;
  for ( i = 0 ; i < 4; i++ )
  {
    for ( j = 0 ; j < 3; j++ )
    {
      if ( flagbox[i][j]->text() == "Hold Previous" )
      {
        if ( hold != flagbox[i][j]->checkState() )
        {
          hold = flagbox[i][j]->checkState();
        };
      }
      if ( flagbox[i][j]->text() == "Hide all" )
      {
        if ( hide_surfaces != flagbox[i][j]->checkState() )
        {
          hide_surfaces = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Reference" )
      {
        if ( drawref != flagbox[i][j]->checkState() )
        {
          drawref = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Approximate" )
      {
        if ( approximate_normals != flagbox[i][j]->checkState() )
        {
          approximate_normals = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Average" )
      {
        if ( average_normals != flagbox[i][j]->checkState() )
        {
          average_normals = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Shiny" )
      {
        if ( shiny != flagbox[i][j]->checkState() )
        {
          shiny = flagbox[i][j]->checkState();
          if ( plotter->NSurf() > 0 )
          {
            int n = plotter->NSurf()-1;
            plotter->surface[n]->shiny = shiny;
            if ( strstr( plotter->surface[n]->name , "Neg" ) != NULL )
            {
              plotter->surface[n-1]->shiny = shiny;
            };
            viewer->Redraw();
          };
        };
      };
    };
  };
}

void DwPlot::SetupPlot()
{
  int table = tabWidget->currentIndex();
  if ( listados[table]->currentItem() == NULL ) return;
  QString st = (listados[table]->currentItem())->text();
  plotter->RequestPlot(st,set_number,orbital_number,perturbation_number,
                       isovalue,true);
}

void DwPlot::BuildPlot()
{
  int table = tabWidget->currentIndex();
  if ( listados[table]->currentItem() == NULL ) 
  {
    gn->ErrorMessage(__FILE__,__LINE__, 
           "You must choose a field first" , 0 );
    return;
  };
  QString st = (listados[table]->currentItem())->text();
  if ( st.isNull() )
  {
  }
  if ( ! hold ) plotter->Clear();
  plotter->RequestPlot(st,set_number,orbital_number,
                       perturbation_number,isovalue,false);
  plotter->RequestPlot(st,set_number,orbital_number,
                       perturbation_number,-isovalue,false);
  viewer->Redraw();
}

void DwPlot::ChangeSetNumber( int n )
{
  set_number = n;
}

void DwPlot::ChangeOrbitalNumber( int n )
{
  orbital_number = n;
}

void DwPlot::ChangePerturbationNumber( int n )
{
  perturbation_number = n;
}

void DwPlot::Adjust()
{
  double l = plotter->adjust_max - plotter->adjust_min;

  cminS->setValue( plotter->adjust_min );
  cmaxS->setValue( plotter->adjust_max );
  cstepS->setValue( 0.05*l );
}

