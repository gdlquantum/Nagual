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
#include <string.h>

#include <iostream>

#include <QtOpenGL/qgl.h>
#include <qcursor.h>
#include <QMessageBox> 
#include <QFileDialog>

#include <GNagual.h>
#include <System.h>
#include <DwPanel.h>
#include <DwElement.h>
#include <Model.h>

using namespace std;
GNagual::GNagual( System* isys, int argc , char* argv[] )
     : QApplication( argc , argv )
{
  sys = isys;

  if ( argc > 2 ) strcpy(input_file,argv[2]);

  panel = new DwPanel(this,sys);

  es = new DwElement(panel,6);

  // Without OpenGL support the program does not work
  if ( !QGLFormat::hasOpenGL() ) 
  {
    qWarning( "This system has no OpenGL support. Exiting." );
    quit();
  };

  panel->show();

  status = 0;
  panel->ChangeModel( QString("Balls & Sticks") );
}

void GNagual::SetCursor( const char *type )
{
  QString typestr = type;
  QCursor cursor;
  if ( typestr == "Normal" )
  {
    cursor = QCursor( Qt::ArrowCursor );
  }
  else if ( typestr == "Wait" )
  {
    cursor = QCursor( Qt::WaitCursor );
  }
  else if ( typestr == "Selection" )
  {
    cursor = QCursor( Qt::PointingHandCursor );
  };
  panel->setCursor( cursor );
} 

void GNagual::ErrorMessage( const char* filename , int linenumber, const char*msg , int el )
{
  QString title = "Attention";
  QString text;
  if ( el )
  {
    text = "Error occurred in ";
  }
  else
  {
    text = "Warning received from "; 
  };
  text += filename;
  text += ", line number ";
  text += QString::number( linenumber );
  text += "\n";
  text += msg;
  text += "\n";

  if ( el )
  {
    (void)QMessageBox::critical( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  }
  else
  {
    (void)QMessageBox::warning( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  };

  if ( el ) closeAllWindows();
}

bool GNagual::GetFileName(char* filename)
{
  QString filter("All Files (*.*)");
  filter +=      "\nOnly    (*.xal)";
  filter +=      "\nOnly    (*.xyz)";
  QString fn = QFileDialog::getOpenFileName(0,
                 QString("Read file"),
                 QDir::currentPath(),
                 filter);
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    sprintf(filename,"%s",(fn.toLatin1()).data());
    return true;
  }
  return false;
}
