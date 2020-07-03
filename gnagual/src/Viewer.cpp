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

#include<iostream>

#include <QtOpenGL/QtOpenGL>
#include <QWheelEvent>

#include <GNagual.h>
#include <Viewer.h>
#include <DwImage.h>
#include <DwPanel.h>
#include <DwGeo.h>
#include <DwPlot.h>
#include <Math.h>
#include <Model.h>
#include <Molecule.h>
#include <Plotter.h>
#include <System.h>
#include <Atom.h>
#include <Element.h>
#include <Units.h>

bool perspective = false;

Viewer::Viewer( GNagual *ig , DwPanel* ipanel )
      : QGLWidget( ipanel )
{
  gn = ig;
  panel = ipanel;
  plotter = new Plotter(gn);
  render = new DwImage( this );

  ResetRotation();

  trans[0] = 0.0; 
  trans[1] = 0.0; 
  trans[2] = -5.0; 
  scale = 1.0;

  bg[0] = 0.0;
  bg[1] = 0.0;
  bg[2] = 0.0;
  bg[3] = 1.0;

  tag_type = "None";
  tag_color[0] = 1.0;
  tag_color[1] = 0.33;
  tag_color[2] = 0.0;
  tag_color[3] = 0.0;

  view_width = width()/10;
  view_height = view_width*((float)height())/(float)width();
  view_near = -100.0;
  view_far = 100.0;
  perspective = false;
  if (perspective) 
  {
    view_near = -10.0;
    view_far = 30.0;
  }

  this->setFocusPolicy(Qt::StrongFocus);

  //rfm user_label.clear();
  selected_atoms.clear();
  monitored_bonds.clear();
  monitored_angles.clear();
  monitored_dihedrals.clear();
}

void Viewer::Project()
{
  makeCurrent();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
}
 
void Viewer::initializeGL( void )
{
  GLfloat whiteDir[4] = {2.0, 2.0, 2.0, 1.0};
  GLfloat lightPos0[4] = {30.0,30.0,30.0,1.0};

  glMaterialf(GL_FRONT, GL_SHININESS, 80.0);
  glLightfv( GL_LIGHT0 , GL_POSITION , lightPos0 );
  glEnable( GL_DEPTH_TEST );
  glEnable( GL_LIGHTING );
  glEnable( GL_LIGHT0 );
  glEnable( GL_COLOR );
  glEnable( GL_COLOR_MATERIAL );
  glShadeModel( GL_SMOOTH );
  glMaterialfv( GL_FRONT, GL_SPECULAR, whiteDir);
  glEnable(GL_NORMALIZE);
  glClearColor( bg[0] , bg[1] , bg[2] , bg[3] );

  resizeGL( width() , height() );
}

void Viewer::resizeGL( int w , int h )
{
  glViewport( 0, 0, w, h );
  view_height = view_width*((float)h/(float)w);
  Project();
}

void Viewer::Redraw()
{
  paintGL();
}

void Viewer::paintGL()
{
  Project();
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );
 
  panel->model->DrawMolecule( -1 , true );

  panel->model->DrawArrows();

  DrawSurfaces();

  DrawStrings();

  glPopMatrix();

  if (context()!=0) swapBuffers();
}

void Viewer::mousePressEvent( QMouseEvent *e )
{
  mouse[0] = e->x();
  mouse[1] = e->y();
  if ( ( gn->status & SFB  ) ||
       ( gn->status & SFA  ) ||
       ( gn->status & SFD  ) ||
       ( gn->status & SFFV ) ||
       ( gn->status & SFRV ) ||
       ( gn->status & SFAV ) ||
       ( gn->status & SFDV ) ||
       ( gn->status & SFPB ) ||       
       ( gn->status & SFPG ) )
  {
    PickAtom();
  };
  if ( gn->status & SFSURF )
  {
    PickSurface();
  }
}

void Viewer::mouseMoveEvent( QMouseEvent *e )
{
  bool rotate = false;
  int sx = e->x() - mouse[0];
  int sy = e->y() - mouse[1];

  if ( abs(sx) > (signed)width()/100 && abs(sx) > abs(sy) ) 
  {
    mouse[0] = e->x();
    Rotate( (double)((180*sx)/(signed)width() % 360), 0.0, 1.0, 0.0); 
    rotate = true;
  };

  if ( abs(sy) > (signed)height()/100 && abs(sx) < abs(sy) ) 
  {
    mouse[1] = e->y();
    Rotate( (double)((180*sy)/(signed)height() % 360), 1.0, 0.0, 0.0); 
    rotate = true;
  };

  if ( rotate )
  {
    Redraw();
  }
}

void Viewer::mouseReleaseEvent( QMouseEvent *e )
{
  mouse[0] = e->x();
  mouse[1] = e->y();
}

void Viewer::keyPressEvent(QKeyEvent *keyEvent)
{
  if (keyEvent->key()==Qt::Key_Up)
  {
    Rotate( -10.0, 1.0, 0.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Down)
  {
    Rotate( 10.0, 1.0, 0.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Left)
  {
    Rotate( -10.0, 0.0, 1.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Right)
  {
    Rotate( 10.0, 0.0, 1.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else
  {
    keyEvent->ignore();
  }
}

void Viewer::Rotate( double angle , double x , double y , double z )
{
  int i,j,k;
  double cosp,sinp,cost,norma;
  double m[3][3];
  double mm[3][3];
  double rm[3][3];
  double u[3];

  u[0] = x;
  u[1] = y;
  u[2] = z;
  norma = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  if ( norma == 0.0 )
  {
    return;
  }
  u[0] = u[0]/norma;
  u[1] = u[1]/norma;
  u[2] = u[2]/norma;
  angle = angle*M_PI/180.0;
  cosp = cos(angle);
  sinp = sin(angle);
  cost = 1.0 - cosp;
  m[0][0] = u[0]*u[0]*cost + cosp;
  m[0][1] = u[0]*u[1]*cost - u[2]*sinp;
  m[0][2] = u[0]*u[2]*cost + u[1]*sinp;
  m[1][0] = u[0]*u[1]*cost + u[2]*sinp;
  m[1][1] = u[1]*u[1]*cost + cosp;
  m[1][2] = u[1]*u[2]*cost - u[0]*sinp;
  m[2][0] = u[0]*u[2]*cost - u[1]*sinp;
  m[2][1] = u[1]*u[2]*cost + u[0]*sinp;
  m[2][2] = u[2]*u[2]*cost + cosp;

  k = 0;
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      rm[j][i] = rot[k];
      k++;
    };
    k++;
  };

  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      mm[i][j] = 0.0;
      for ( k = 0 ; k < 3 ; k++ )
      {
        mm[i][j] += m[i][k]*rm[k][j];
      };
    };
  };

  k = 0;
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      rot[k] = mm[j][i];
      k++;
    };
    k++;
  };

}
void Viewer::ResetRotation(void)
{
  for ( int i = 0 ; i < 16 ; i++ )
    rot[i] = 0.0;
  rot[0]  = 1.0;
  rot[5]  = 1.0;
  rot[10] = 1.0;
  rot[15] = 1.0;
}

void Viewer::ZoomIn()
{
  scale *= 1.1;
  Redraw();
}

void Viewer::ZoomOut()
{
  scale /= 1.1;
  Redraw();
}

void Viewer::wheelEvent(QWheelEvent *event)
{
    if (event->delta()>0){
        ZoomIn();
    }else{
        if (event->delta()<0){
            ZoomOut();
        }
    }
}

void Viewer::PickSurface( void )
{
  GLfloat picksize=3.0;
  GLuint buffer[512];
  GLint viewport[4];
  GLint hits;

  glGetIntegerv( GL_VIEWPORT , viewport );
  glSelectBuffer( 512 , buffer );
  glInitNames();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPickMatrix( (GLfloat)mouse[0] ,
                 (GLfloat)(viewport[3]-mouse[1]) ,
                  picksize , picksize , viewport );
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  glRenderMode( GL_SELECT );
  glMatrixMode( GL_MODELVIEW );
  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );

  for ( int i = 0 ; i < plotter->NSurf() ; i++ )
  {
    glPushMatrix();
    glPushName( i );
    plotter->DoSurface( i );
    glPopName();
    glPopMatrix();
  };

  glPopMatrix();
  hits = glRenderMode( GL_RENDER );
  Project();
  if ( hits > 0 )
  {
    SaveSelection( (int) buffer[3] );
  }; 
}

int Viewer::SelectedAtom( size_t n )
{
  for ( size_t i = 0 ; i < selected_atoms.size() ; i++ )
  {
    if ( n == selected_atoms[i] )
    {
      return i; 
    };
  };
  return -1; 
}

void Viewer::SaveSelection( int n )
{
  bool changed;

  changed = false;

  if ( gn->status & SFSURF )
  {
    gn->beep();
    gn->status &= ~SFSURF;
    QColor color; 
    color = QColorDialog::getColor( color , this );
    double col[4];
    plotter->GetColor( n , col );
    col[0] = color.red()/255.0;
    col[1] = color.green()/255.0;
    col[2] = color.blue()/255.0;
    plotter->SetColor( n , col );
    gn->SetCursor( "Normal" );
    Redraw();
  }
  else if ( gn->status & SFRV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      gn->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 2 )
      {
        monitored_bonds.push_back( selected_atoms );
        selected_atoms.clear();
        gn->status &= ~SFRV;
        panel->Report("Atom %d selected, selection done",n + 1);  
        gn->SetCursor( "Normal" );
      }
      else
      {
        panel->Report("Atom %d selected, select another",n + 1);  
      };
      Redraw();
    };
  }

 else if ( gn->status & SFPG )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      gn->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 3 )
      {
        int ref[3];
        ref[0] = selected_atoms[0];
        ref[1] = selected_atoms[1];
        ref[2] = selected_atoms[2];
        gn->status &= ~SFPG;
        if (panel->geoeditor->SGDrv(ref))
        {
          panel->Report("Group building failed");
        }
        else
        {
          panel->geoeditor->Update();
          changed = true;
        }
      }
      else
      {
        Redraw();
      };
    };
  }
 else if ( gn->status & SFPB )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      gn->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 2 )
      {
        gn->status &= ~SFPB;
        changed = true;
      }
      else
      {
        panel->Report("Atom %d selected, select another",n + 1);
      };
      Redraw();
    };
  }
  else if ( gn->status & SFAV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      gn->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 3 )
      {
        monitored_angles.push_back( selected_atoms );
        selected_atoms.clear();
        gn->status &= ~SFAV;
        panel->Report("Atom %d selected, selection done",n + 1);  
        gn->SetCursor( "Normal" );
      }
      else
      {
        panel->Report("Atom %d selected",n + 1);  
      };
      Redraw();
    };
  }
  else if ( gn->status & SFDV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      gn->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 4 )
      {
        monitored_dihedrals.push_back( selected_atoms );
        selected_atoms.clear();
        gn->status &= ~SFDV;
        panel->Report("Atom %d selected, selection done",n + 1);  
        gn->SetCursor( "Normal" );
      }
      else
      {
        panel->Report("Atom %d selected",n + 1);  
      };
      Redraw();
    };
  } 
  else if ( gn->status & SFB )
  {
    gn->beep();
    panel->Report( "Atom %d selected",n+1);
    gn->status &= ~SFB;
    if ( panel->geoeditor->mol->Natom() == 2 )
    {
      panel->geoeditor->AddAtom( n + 1, 2 - n, 1);
      changed = true;
    }
    else
    {
      gn->status |= SFA;
      selected_atoms.push_back( n );
      Redraw();
      panel->Report("Atom %d selected",n + 1);
    };
  }
  else if ( gn->status & SFA )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      selected_atoms.push_back( n );
      gn->beep();
      panel->Report( "Atom %d selected",n+1);
      gn->status &= ~SFA;
      if ( panel->geoeditor->mol->Natom() == 3 )
      {
        panel->geoeditor->AddAtom( selected_atoms[0] + 1, n + 1,
                                   4 - selected_atoms[0] - n);
        changed = true;
      }
      else
      {
        gn->status |= SFD;
        Redraw();
        panel->Report( "Atom %d selected",n + 1);
      };
    };
  }
  else if ( gn->status & SFD )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      selected_atoms.push_back( n );
      gn->beep();
      panel->Report("Atom %d selected",n+1);
      gn->status &= ~SFD;
      panel->geoeditor->AddAtom(selected_atoms[0] + 1,
                                selected_atoms[1] + 1, n + 1);
      changed = true;
    };
  };

  if ( changed )
  {
    selected_atoms.clear();
    gn->SetCursor( "Normal" );
  };
}

void Viewer::DrawSurfaces()
{
  for ( int i = 0 ; i < plotter->NSurf(); i++ )
  {
    plotter->DoSurface( i );
  };
};

void Viewer::Write(char* text, Vector p, size_t fontsize)
{
  QFont font( "Helvetica" , fontsize );
  glColor3f( tag_color[0] , tag_color[1] , tag_color[2] );
  glDisable( GL_LIGHTING );
  glEnable( GL_LINE_STIPPLE );
  glLineStipple( 1 , 0xf0f0 );
  renderText( p[0] , p[1] , p[2] , QString(text), font );
  glDisable( GL_LINE_STIPPLE );
  glEnable( GL_LIGHTING );
}

void Viewer::ChangeBackgroundColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  makeCurrent();
  bg[0] = color.red()/255.0;
  bg[1] = color.green()/255.0;
  bg[2] = color.blue()/255.0;
  glClearColor( bg[0] , bg[1] , bg[2] , bg[3] );
  Redraw();
}

void Viewer::ChangeTags( const QString &new_tag_type )
{
  tag_type = new_tag_type;
  Redraw();
}

void Viewer::ChangeTagColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  tag_color[0] = color.red()/255.0;
  tag_color[1] = color.green()/255.0;
  tag_color[2] = color.blue()/255.0;
  Redraw();
}


void Viewer::PickAtom( void )
{
  GLfloat picksize=3.0;
  GLuint buffer[512];
  GLint viewport[4];
  GLint hits;
  int iatom;

  glGetIntegerv( GL_VIEWPORT , viewport );
  glSelectBuffer( 512 , buffer );
  glInitNames();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPickMatrix( (GLfloat)mouse[0] ,
                 (GLfloat)(viewport[3]-mouse[1]) ,
                 picksize , picksize , viewport );
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  glRenderMode( GL_SELECT );
  glMatrixMode( GL_MODELVIEW );
  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );
  for ( iatom = 0 ; iatom < panel->geoeditor->mol->Natom() ; iatom++ )
  {
    glPushMatrix();
    glPushName( iatom );
    panel->model->DrawMolecule( iatom , false );
    glPopName();
    glPopMatrix();
  }

  glPopMatrix();
  hits = glRenderMode( GL_RENDER );
  Project();
  if ( hits > 0 )
  {
    SaveSelection( (int) buffer[3] );
  }; 
}

void Viewer::DrawStrings()
{
  static int angle_step = 3;
  size_t i,j,iatom,jatom,katom;
  size_t fontsize=32;
  char word[GN_MAX_STR_SIZE];
  double val,ss;
  Vector va,vb,vc,vd,vshift;
  Vector rva,rvb,rvc,rvd;
  Vector auxa,auxb,auxc;

  if ( bg[0] == 1.0 &&  bg[1] == 1.0 && bg[2] == 1.0 )
   glColor3f( 0.0 , 0.0 , 0.0 );
  else glColor3f( 1.0 , 1.0 , 1.0 );

  vshift = Vector( rot[2], rot[6] , rot[10] );

  glDisable( GL_LIGHTING );
  glEnable( GL_LINE_STIPPLE );
  glLineStipple( 1 , 0xf0f0 );

  //rfm strcpy( word , (image_render->label->text().toAscii()).data() );
  fontsize = 24;
  //rfm renderText( 10 , 50 , QString(word), QFont( "Helvetica" , fontsize ) );
  ss = 5.0*panel->model->crs*scale;
  fontsize = size_t(GN_MAX(4.0,18*scale));
/*
  for ( i = 0 ; i < monitored_bonds.size() ; i++ )
  {
    iatom = monitored_bonds[i][0];
    jatom = monitored_bonds[i][1];
    va = Vector( &COORD[1][iatom][0] );
    vb = Vector( &COORD[1][jatom][0]);
    glBegin( GL_LINES );
      glVertex3d( va[0] , va[1] , va[2] );
      glVertex3d( vb[0] , vb[1] , vb[2] );
    glEnd();
    vc = va - vb;
    val = BohrToAngstrom(vc.Norm()); 
    va += vb;
    va *= 0.5;
    vb = vshift;
    vb *= ss;
    va += vb;
    sprintf(word,"%3.3f",val);
    Write(word, va, fontsize );
  };

  for ( i = 0 ; i < monitored_angles.size() ; i++ )
  {
    iatom = monitored_angles[i][0];
    jatom = monitored_angles[i][1];
    katom = monitored_angles[i][2];
    va = Vector(&COORD[1][iatom][0]);
    vb = Vector(&COORD[1][jatom][0]);
    vc = Vector(&COORD[1][katom][0]);

    rva = va - vb;
    rvc = vc - vb;
    rva ^= 1.0;
    rvc ^= 1.0;
    rvb = rva > rvc;
    val = 180.0/M_PI*acos( rva.Dot( rvc ) );
    sprintf(word,"%3.1f",val);
    rva ^= 1.0;
    rva.Rotate( rvb , val/2.0 );
    rvc = vb + rva;
    rva.Rotate( rvb , -val/2.0 );
    Write(word, rvc, fontsize );

    glBegin( GL_LINES );
    glVertex3d( va[0] , va[1], va[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    glEnd();
    glBegin( GL_POINTS );
    for ( j = 0 ; j < (size_t) val ; j += angle_step )
    {
      rva.Rotate( rvb , (double)angle_step );
      rvc = vb + rva;
      glVertex3d( rvc[0] , rvc[1] , rvc[2] );
    };
    glEnd();
  };

  for ( i = 0 ; i < monitored_dihedrals.size() ; i++ )
  {
    iatom = monitored_dihedrals[i][0];
    va = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][1];
    vb = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][2];
    vc = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][3];
    vd = Vector(&COORD[1][iatom][0]);

    rva = va - vb;
    rvc = vc - vb;
    rvd = vd - vc;
    rva ^= 1.0;
    rvc ^= 1.0;
    rvd ^= 1.0;
    auxa = rvc;
    auxa *= -rva.Dot( rvc );
    auxa += rva;
    auxa ^= 1.0;
    auxb = rvc;
    auxb *= -rvd.Dot( rvc );
    auxb += rvd;
    auxb ^= 1.0;
    auxc = auxa;
    auxc.Rotate( rvc , -90.0 );
    val = auxa.Dot( auxb );
    if ( val < 0.0 )  
    {
      if ( auxb.Dot( auxc ) < 0.0 )
      {
        val = 360.0 - 180.0/M_PI*acos( val );
      }
      else
      {
        val =  180.0/M_PI*acos( val );
      }
    }
    else
    {
      if ( auxb.Dot( auxc ) < 0.0 )
      {
        val = 360.0 - 180.0/M_PI*acos( val );
      }
      else
      {
        val = 180.0/M_PI*acos( val );
      }
    };
    sprintf(word,"%3.1f",val);
    auxa = rvc;
    auxa *= -rva.Dot( rvc );
    auxa += rva;
    auxa.Rotate( rvc , -val/2.0 );
    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    auxa.Rotate( rvc , val/2.0 );
    Write(word, auxb,  fontsize);

    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    glBegin( GL_LINES );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( va[0] , va[1], va[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    auxa.Rotate( rvc , -val );
    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    auxa.Rotate( rvc , val );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vd[0] , vd[1], vd[2] );
    glEnd();

    glBegin( GL_POINTS );
    for ( j = 0 ; j < (size_t) val ; j += angle_step )
    {
      auxa.Rotate( rvc , -(double)angle_step );
      auxb = vb + vc; 
      auxb *= 0.5; 
      auxb += auxa;
      glVertex3d( auxb[0], auxb[1], auxb[2] );
    };
    glEnd();
  };

  for ( i = 0 ; i < user_label.size() ; i++ )
  {
    renderText( user_label[i].x, 
                user_label[i].y, 
                QString( user_label[i].text ), QFont( "Helvetica" , 32 )); 
  }
*/

  // Atom labels are orange (Roberto Flores-Moreno, Feb 2008)
  glColor3f( tag_color[0] , tag_color[1] , tag_color[2] );
  if ( ! (tag_type == "None" ) )
  {
    Atom *atom;
    ss = 10.0*panel->model->srs*scale;
    vb = vshift;
    vb *= ss;
    for ( iatom = 0 ; iatom < panel->geoeditor->mol->Natom(); iatom++ )
    {
      atom = panel->geoeditor->mol->atom[iatom];
      va = Vector(atom->x,atom->y,atom->z);
      va += vb;
      if ( tag_type == "Symbols" )
      {
        sprintf(word,"%s",ELEMENT_SYMBOL[atom->atomic_number]);
      }
      else if ( tag_type == "Numbers" )
      {
        sprintf(word,"(%d)",iatom+1);
      }
      else if ( tag_type == "Sym+Num" )
      {
        sprintf(word,"%s(%d)",ELEMENT_SYMBOL[atom->atomic_number],iatom+1);
      }
      else if ( tag_type == "Charges" )
      {  
        sprintf(word,"%.2f",atom->charge);
      }
/*
      else if ( tag_type == "Fukui HOMO" )
      {  
        sprintf(word,"%.3f",ATOMFUKUI_H[iatom]);
      }
      else if ( tag_type == "Fukui LUMO" )
      {  
        sprintf(word,"%.3f",ATOMFUKUI_L[iatom]);
      }
      else if ( tag_type == "Fukui Average" )
      {  
        sprintf(word,"%.3f",(ATOMFUKUI_H[iatom]+ATOMFUKUI_L[iatom])/2.0);
      }
      else if ( tag_type == "Fukui Difference" )
      {  
        sprintf(word,"%.2f",(ATOMFUKUI_L[iatom]-ATOMFUKUI_H[iatom]));
      }
*/
      else
      {
        break;
      }
      Write( word,va, fontsize );
    }
  };

  glEnable( GL_LIGHTING );
  glDisable( GL_LINE_STIPPLE );
}



