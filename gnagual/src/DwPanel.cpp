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
#include <stdarg.h>
#include <pthread.h>

#include <iostream>
#include <fstream>

#include <QtCore>
#include <QtGui>
#include <QMainWindow>

#include <DwPanel.h>   
#include <GNagual.h>   
#include <Viewer.h>   
#include <DwImage.h>   
#include <DwGeo.h>   
#include <DwPlot.h>   
#include <System.h>   
#include <Molecule.h>   
#include <Model.h>   
#include <DwElement.h>
#include <Element.h>

using namespace std;

static const char *READ_FORMATS[] = {"NGL","XYZ",NULL};

DwPanel::DwPanel( GNagual *ig, System* isys )
     : QWidget( 0 )
{
  gn = ig;
  sys = isys;

  tabWidget = new QTabWidget;

  viewer = new Viewer( gn , this );
  viewer->hide();
  geoeditor = new DwGeo( this );
  geoeditor->hide();
  model = new Model( this );
  winplot = new DwPlot( gn , viewer );


  statusbar = new QStatusBar( this );

 int i;


  i = 0;
  while ( READ_FORMATS[i] != NULL )
  {
    RFMTA.push_back( new QAction( tr( READ_FORMATS[i] ), this ) );
    RFMTA[i]->setStatusTip(tr("Open file in %1 format")
                           .arg(READ_FORMATS[i]));
    RFMTA[i]->setCheckable( true );
    connect( RFMTA[i], SIGNAL(triggered()), this, SLOT( ReadFile() ));
    i++;
  };

  imgAct = new QAction(tr("&Image"), this);
  imgAct->setStatusTip(tr("Write image file"));
  connect(imgAct, SIGNAL(triggered()), this, SLOT(WriteImage()));

  quitAct = new QAction(tr("&Quit"), this);
  quitAct->setShortcut(tr("Ctrl+Q"));
  quitAct->setStatusTip(tr("Finish the application"));
  connect(quitAct, SIGNAL(triggered()), gn, SLOT(closeAllWindows()));

  QMenuBar *mainMenu = new QMenuBar;

  fileMenu = new QMenu(tr("&File"),this);

  readMenu = new QMenu(tr("&Read"),this);
  for( size_t i = 0 ; i < RFMTA.size() ; i++ )
    readMenu->addAction( RFMTA[i] );
  fileMenu->addMenu(readMenu);

  fileMenu->addAction( imgAct );

  fileMenu->addSeparator();
  fileMenu->addAction( quitAct );

  mainMenu->addMenu(fileMenu);

  QGridLayout *mainLayout = new QGridLayout;

  inButton = new QToolButton( this );
  inButton->setText( "Zoom In" );
  connect( inButton,SIGNAL(clicked()),viewer,SLOT(ZoomIn()));
  
  outButton = new QToolButton( this );
  outButton->setText( "Zoom Out" );
  connect( outButton,SIGNAL(clicked()),viewer,SLOT(ZoomOut()));

  giraButton = new QToolButton( this );
  giraButton->setText( "Start Rotation" );
  connect( giraButton,SIGNAL(clicked()),this,SLOT(SetSpin()));

  modelComboBox = new QComboBox;
  modelComboBox->addItem( "Model" );
  modelComboBox->addItem( QString("None"));
  modelComboBox->addItem( QString("Balls & Sticks"));
  modelComboBox->addItem( QString("Balls"));
  modelComboBox->addItem( QString("Sticks"));
  modelComboBox->addItem( QString("Wireframe"));
  modelComboBox->addItem( QString("Van der Waals"));
  modelComboBox->addItem( QString("Custom"));
  connect(modelComboBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeModel(const QString &)));

  QComboBox *tagComboBox = new QComboBox;
  tagComboBox->addItem( "None");
  tagComboBox->addItem( "Symbols");
  tagComboBox->addItem( "Numbers");
  tagComboBox->addItem( "Sym+Num");
  tagComboBox->addItem( "Charges");
  connect(tagComboBox, SIGNAL(activated(const QString &)),
          viewer, SLOT(ChangeTags(const QString &)));
  QLabel* taglab = new QLabel(tr("Labels:"));
  taglab->setBuddy(tagComboBox);


  arrowComboBox = new QComboBox;
  arrowComboBox->addItem( "None");
  arrowComboBox->addItem( "Forces");
  connect(arrowComboBox, SIGNAL(activated(const QString &)),
          model, SLOT(ChangeArrows(const QString &)));
  QLabel* arrowlab = new QLabel(tr("Arrows:"));
  arrowlab->setBuddy(tagComboBox);

  styleComboBox = new QComboBox;
  styleComboBox->addItems(QStyleFactory::keys());
  connect(styleComboBox, SIGNAL(activated(const QString &)),
          this, SLOT(MakeStyle(const QString &)));


  QToolButton *writeButton =  new QToolButton( this );
  writeButton->setText( "Write" );
  connect(writeButton,SIGNAL(clicked()),this,SLOT(WriteFile()));

  QWidget *visualizer = new QWidget(0);
  QVBoxLayout *visuallayout = new QVBoxLayout;

  visuallayout->addWidget( modelComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( inButton,1,Qt::AlignCenter);
  visuallayout->addWidget( outButton,1,Qt::AlignCenter);
  visuallayout->addWidget( giraButton,1,Qt::AlignCenter);
  visuallayout->addWidget( styleComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( taglab,1,Qt::AlignCenter);
  visuallayout->addWidget( tagComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( arrowlab,1,Qt::AlignCenter);
  visuallayout->addWidget( arrowComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( writeButton,1,Qt::AlignCenter);
  visuallayout->addWidget( statusbar,1,Qt::AlignCenter);

  styleComboBox->hide();

  visualizer->setLayout(visuallayout);

  mainLayout->setMenuBar(mainMenu);
  mainLayout->addWidget( viewer        ,      0 , 0 , 22 , 20);
  mainLayout->addWidget( tabWidget     ,      0 ,20 ,  0 ,  2);
  mainLayout->addWidget( statusbar     ,      22, 0 ,  1 , 20);

  setLayout( mainLayout );

  tabWidget->clear();
  tabWidget->addTab( visualizer , QString("Control") );
  tabWidget->addTab( winplot , QString("Plotter") );
  tabWidget->addTab( geoeditor , QString("Geometry") );

  timer = new QTimer(this);

  setWindowTitle( "GNagual 0.7" );
  move( 0 , 0  );
  resize(1500 , 1000);
  viewer->show();
}


void DwPanel::Report( const char *fmt, ... )
{
  va_list args;   
  char    s[256];

  va_start( args , fmt );
  vsprintf(s , fmt , args );
  va_end( args );
 
  if (statusbar) statusbar->showMessage( QString( s ) );
}

void DwPanel::About()
{
  string str = "GNagual software\n";
  str += "Authors: \n";
  str += "R. Flores-Moreno      roberto.floresmoreno.qt@gmail.com\n";
  QMessageBox::about( this, "About GNagual",str.c_str());
}

void DwPanel::MakeStyle( const QString &style )
{
  QApplication::setStyle(QStyleFactory::create(style));
}

void DwPanel::ReadFile()
{
  char fileformat[256];
  size_t i;

  for ( i = 0 ; i < RFMTA.size() ; i++ )
  {
    if ( RFMTA[i]->isChecked() ) 
    {
      RFMTA[i]->setChecked( false );
      sprintf(fileformat,"%s",((RFMTA[i]->text()).toLatin1()).data());
      break;
    };
  };

  QString filter("All Files (*.*)");
  filter +=      "\nOnly    (*.xal)";
  filter +=      "\nOnly    (*.xyz)";
  QString fn = QFileDialog::getOpenFileName(this,
                 QString("Read file"),
                 QDir::currentPath(),
                 filter);
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[256];

    sprintf(filename,"%s",(fn.toLatin1()).data());
    Report( "Reading file %s ...", filename );
    string tit = "Green";
    tit += ":";
    tit += filename;
    setWindowTitle( tit.c_str() );
    DoReadFile(filename,fileformat);
    Report( "Reading file %s ... DONE", filename );
  };
}

void DwPanel::DoReadFile(char* filename, char* fmt)
{
  sys->Read(filename,fmt);
  geoeditor->Update();
  model->ChangeKind("Balls & Sticks");
  viewer->Redraw();
}

void DwPanel::WriteFile()
{
  char sffmt[GN_MAX_STR_SIZE];
  char fileformat[GN_MAX_STR_SIZE];
  QString ffmt;
  size_t i;

  strcpy(fileformat,"XYZ");

  QString initialPath = QDir::currentPath() + "/untitled.";
  initialPath = initialPath + ffmt.toLower();
  QString fn = QFileDialog::getSaveFileName( this,
               tr( "Write %1 file").arg(ffmt.toUpper()),
               initialPath, tr("All Files (*.%1)").arg(ffmt.toLower()) );
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[GN_MAX_STR_SIZE];

    sprintf(filename,"%s",(fn.toLatin1()).data());
    Report( "Writing file %s ...", filename );
    DoWriteFile(filename , fileformat );
    Report( "Writing file %s ... DONE", filename );
  };
}

void DwPanel::DoWriteFile(char* filename, char* fmt)
{
  if ( strcasecmp( fmt , "XYZ" ) == 0 )
  {
    geoeditor->mol->WriteXYZ(filename,true);
  }
  else if ( strcasecmp( fmt , "ZMT" ) == 0 )
  {
    geoeditor->mol->WriteZMatrix(filename);
  }
}

void DwPanel::WriteImage()
{ 
  viewer->render->show(); 
}

void DwPanel::ChangeAtomColor()
{
  gn->es->exec(); 
  int sel = gn->es->GetSelectedElement();
  QColor color; 
  color = QColorDialog::getColor( color , this );
  ELEMENT_COLOR[sel][0] = color.red()/255.0;
  ELEMENT_COLOR[sel][1] = color.green()/255.0;
  ELEMENT_COLOR[sel][2] = color.blue()/255.0;
  viewer->Redraw();
}

void DwPanel::ChangeSurfaceColor()
{
  gn->SetCursor( "Selection" );
  gn->status |= SFSURF;
  Report("Select a surface");
}


void DwPanel::SetSpin(void)
{
  if (giraButton->text()=="Start Rotation")
  {
    giraButton->setText("Stop Rotation");
    connect(timer,SIGNAL(timeout()),
            this,SLOT(AdvanceSpin()));
    AdvanceSpin();
  }
  else if (giraButton->text()=="Stop Rotation")
  {
    giraButton->setText("Start Rotation");
    disconnect(timer,SIGNAL(timeout()),
               this,SLOT(AdvanceSpin()));
  }
}

void DwPanel::AdvanceSpin(void)
{
  viewer->Rotate( 5.0, 0.0, 1.0, 0.0); 
  viewer->Redraw();
  timer->start(40);
}

void DwPanel::ChangeModel( const QString &newmodel )
{
  model->ChangeKind(newmodel);
}


