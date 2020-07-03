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
#include <QtGui>

#include <QtOpenGL/QGLWidget>
#include <DwImage.h>
#include <DwPanel.h>
#include <Viewer.h>

DwImage::DwImage( Viewer *iviewer )
       : QWidget( 0 ) 
{
        viewer = iviewer;
        QVBoxLayout *mainLayout = new QVBoxLayout;
//
        QGroupBox *labelGroupBox = new QGroupBox(tr("Label"));
        QHBoxLayout *labellayout = new QHBoxLayout;
        label = new QLineEdit( labelGroupBox ); 
        label->setText( QString( "" ) ); 
       // connect( label , SIGNAL( returnPressed() ),
       //          this , SLOT( RedrawGLWidget() ));
        labellayout->addWidget( label );
        labelGroupBox->setLayout( labellayout );
        mainLayout->addWidget( labelGroupBox );
//
        QGroupBox *horizontalGroupBox = new QGroupBox(tr("Axis"));
        QHBoxLayout *axislayout = new QHBoxLayout;
        for (int i = 0; i < 3; i++) 
        {
          axised[i] = new QLineEdit( horizontalGroupBox ); 
          axised[i]->setText( QString::number( 1.0 ) ); 
          axislayout->addWidget( axised[i] );
        };
        horizontalGroupBox->setLayout( axislayout );
        mainLayout->addWidget( horizontalGroupBox );
//
        QHBoxLayout *descLayout = new QHBoxLayout;
        formatCombo = new QComboBox;
        formatCombo->addItem( QString( "PNG" ));
        formatCombo->addItem( QString( "BMP" ));
        descLayout->addWidget( formatCombo );
        connect( formatCombo, SIGNAL(activated(const QString &)),
                 this, SLOT( ChangeFormat(const QString &)));
//
        QPushButton *writeB = new QPushButton( tr("Write") );
        connect( writeB, SIGNAL(clicked()), this, SLOT(Write()));
        descLayout->addWidget( writeB );
//
        QPushButton *doneB = new QPushButton( tr("Done") );
        connect( doneB, SIGNAL(clicked()), this, SLOT(close()));
        descLayout->addWidget( doneB );
//
        mainLayout->addLayout( descLayout );
        setLayout( mainLayout );
//
        setWindowTitle( "Image Renderer" );

}

DwImage::~DwImage()
{
}

bool DwImage::GetFileName( char* format )
{
    //    QString initialPath = QDir::currentPath() + "/untitled." + format;
  QString initialPath = QDir::currentPath();
  QString fn = QFileDialog::getSaveFileName(this, 
   "Save Image to File", initialPath,format);
  if ( fn.isEmpty() || fn.isNull() ) 
  {
    return false;
  }
  else
  {
    sprintf(filename,"%s",(fn.toAscii()).data());
    return true;
  }
}

void DwImage::ChangeFormat( const QString& newfmt )
{
  axised[0]->setEnabled( false );
  axised[1]->setEnabled( false );
  axised[2]->setEnabled( false );
}

void DwImage::Write( )
{
  QString text = formatCombo->currentText();
  char format[64];
  sprintf( format , "%s", ((text.toLower()).toAscii()).data() );
  if ( GetFileName( format ) )
  {
    viewer->panel->Report( "Writing %s image file %s ...",
                            format ,  filename );
    QPixmap pixmap = viewer->renderPixmap(viewer->width(), 
                                          viewer->height());
    pixmap.save( filename , format , 0 );
    viewer->panel->Report( "Writing %s image file %s ... DONE",
                           format , filename);
  };
}

