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
#include <GL/glu.h>

#include <Model.h>
#include <DwPanel.h>
#include <Molecule.h>
#include <System.h>
#include <Atom.h>
#include <Viewer.h>

Model::Model(DwPanel* p)
{
  panel = p;
  mol = panel->sys->mol;

  srs = 1.0;
  crs = 1.0;
  drawmol = true;

  arrow_type = "None";
  arrow_color[0] = 0.0;
  arrow_color[1] = 0.0;
  arrow_color[2] = 0.0;
  arrow_color[3] = 0.0;

  q = gluNewQuadric();
};

void Model::ChangeKind( const QString &newmodel )
{
  if ( newmodel == "Model" ) return;

  if ( newmodel == "None" ) drawmol = false;
  else drawmol = true;

  if ( newmodel == "Wireframe" )
  {
    srs = 0.0;
    crs = 0.0;
  }
  else if ( newmodel == "Balls" )
  {
    srs = 0.3;
    crs = 0.0;
  }
  else if ( newmodel == "Balls & Sticks" )
  {
    srs = 0.2;
    crs = 0.12;
  }
  else if ( newmodel == "Sticks" )
  {
    srs = -0.1;
    crs = -0.1;
  }
  else if ( newmodel == "Van der Waals" )
  {
    srs = -0.0001;
    crs = 0.0;
  };

  panel->viewer->Redraw();
}


void Model::DrawSphere( double* p , double radi , int n )
{
  glTranslatef ( p[0] , p[1] , p[2] );
  gluSphere( q , radi , n, n/2 );
  glTranslatef (-p[0] ,-p[1] ,-p[2] );
}

void Model::DrawCylinder(double *start, double *end, 
                          double radi , int n )
{
  GLdouble v[3];
  GLdouble dval, ca, sa, cb, sb;
  GLdouble tmat[4][4];
  GLfloat height; 

  v[0] = (GLdouble)(end[0]-start[0]); 
  v[1] = (GLdouble)(end[1]-start[1]); 
  v[2] = (GLdouble)(end[2]-start[2]); 
  height = (GLfloat) sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  dval = sqrt(v[0]*v[0] + v[1]*v[1]); 
  if (dval != 0.0) 
  {
    ca = v[1]/dval;
    sa = v[0]/dval;
    v[1] = v[0]*sa + v[1]*ca;
  } 
  else
  {
    ca = 1.0;
    sa = 0.0;
  };

  dval = sqrt(v[1]*v[1] + v[2]*v[2]); 
  if (dval != 0.0) 
  {
    cb = v[2]/dval;
    sb = v[1]/dval;
  }
  else 
  {
    cb = 1.0;
    sb = 0.0;
  };

  tmat[0][0] = ca;   
  tmat[1][0] = sa*cb; 
  tmat[2][0] = sa*sb; 
  tmat[3][0] = (GLdouble)start[0];
  tmat[0][1] = (-sa); 
  tmat[1][1] = ca*cb; 
  tmat[2][1] = ca*sb; 
  tmat[3][1] = (GLdouble)start[1];
  tmat[0][2] = 0.0;   
  tmat[1][2] = (-sb); 
  tmat[2][2] = cb;    
  tmat[3][2] = (GLdouble)start[2];
  tmat[0][3] = 0.0;   
  tmat[1][3] = 0.0;   
  tmat[2][3] = 0.0;   
  tmat[3][3] = 1.0;

  glPushMatrix();
  glMultMatrixd((GLdouble*)tmat);
  gluCylinder( q, (GLfloat)radi, (GLfloat)radi, height, n, 1);
  glPopMatrix();
}

void Model::DrawAtom( Atom* a, double radi , int n )
{
  double p[3];

  p[0] = a->x;
  p[1] = a->y;
  p[2] = a->z;

  glColor3dv( a->GetColor() );

  DrawSphere( p , radi , n );
}

void Model::DrawBond(Atom* a , Atom* b , double radi , int n )
{
  double start[3],end[3],mid[3];

  start[0] = a->x;
  start[1] = a->y;
  start[2] = a->z;
  end[0] = b->x;
  end[1] = b->y;
  end[2] = b->z;
  mid[0] = (start[0]+end[0])/2.0;
  mid[1] = (start[1]+end[1])/2.0;
  mid[2] = (start[2]+end[2])/2.0;

  gluQuadricDrawStyle( q , GLU_FILL );
  glColor3dv( a->GetColor() );
  DrawCylinder( start , mid , radi , n );

  glColor3dv( b->GetColor() );
  DrawCylinder( mid, end , radi , n );
}


void Model::DrawMolecule( int id, bool bonds)
{
  int iatom,jatom;
  int ns;
  double dist,radi;
  Atom *a,*b;
  double vec1[4],vec2[4];

  if (!drawmol) goto ARROW;

  ns = 50;
  glEnable( GL_LIGHTING );
  if ( srs == 0.0 && crs == 0.0 )
  {
    glDisable( GL_LIGHTING );
    glLineWidth( 2 );
    glBegin( GL_LINES );
  };

  int ll,ul;
  if (id<0)
  {
    ll = 0;
    ul = mol->Natom()-1; 
  }
  else
  {
    ll = id;
    ul = id;
  }
  for ( iatom = ll ; iatom <= ul ; iatom++ )
  {
    a = mol->GetAtom(iatom);
    if ( srs != 0.0 )
    {
      if ( srs < 0.0 )
      {
        if ( srs == -0.0001 )
        {
          radi = a->GetVDWRadius();
        }
        else
        {
          radi = -srs;
        }
      }
      else
      {
        radi = srs*a->GetVDWRadius();
      };
      if ( a->atomic_number == 0 ) 
      {
        radi = 0.10;
      };
      DrawAtom( a , radi , ns );
    };

    if (crs != 0.0 || srs == 0.0 ) 
    {
      for ( jatom = iatom-1 ; jatom >= 0 ; jatom-- )
      {
        b = mol->GetAtom(jatom);
        dist = b->Distance(a);
        if ( dist < 1.3*(a->GetCovalentRadius() + b->GetCovalentRadius())&&
          a->atomic_number != 0 && b->atomic_number != 0 )
        {
          if ( srs == 0.0 && crs == 0.0 )
          {
            glColor3dv( a->GetColor() );
            glVertex3d( a->x, a->y, a->z ); 
            glColor3dv( b->GetColor() );
            glVertex3d( b->x, b->y, b->z ); 
          }
          else
          {
            if ( crs < 0.0 )
            {
              radi = - crs;
            }
            else
            {
              radi = a->GetCovalentRadius();
              if ( radi > b->GetCovalentRadius() )
              {
                radi = b->GetCovalentRadius();
              };
              radi *= crs;
            };
            DrawBond( a , b , radi , ns/2 ); 
          };
        };
      };
    };
  };
ARROW:

  if ( srs == 0.0 && crs == 0.0 )
  {
    glEnable( GL_LIGHTING );
    glLineWidth( 1 );
    glEnd();
  };
} 

void Model::DrawCone(double *start, double *end, 
                          double radi , int n )
{
  GLdouble v[3];
  GLdouble dval, ca, sa, cb, sb;
  GLdouble tmat[4][4];
  GLfloat height; 

  v[0] = (GLdouble)(end[0]-start[0]); 
  v[1] = (GLdouble)(end[1]-start[1]); 
  v[2] = (GLdouble)(end[2]-start[2]); 
  height = (GLfloat) sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  dval = sqrt(v[0]*v[0] + v[1]*v[1]); 
  if (dval != 0.0) 
  {
    ca = v[1]/dval;
    sa = v[0]/dval;
    v[1] = v[0]*sa + v[1]*ca;
  } 
  else
  {
    ca = 1.0;
    sa = 0.0;
  };

  dval = sqrt(v[1]*v[1] + v[2]*v[2]); 
  if (dval != 0.0) 
  {
    cb = v[2]/dval;
    sb = v[1]/dval;
  }
  else 
  {
    cb = 1.0;
    sb = 0.0;
  };

  tmat[0][0] = ca;   
  tmat[1][0] = sa*cb; 
  tmat[2][0] = sa*sb; 
  tmat[3][0] = (GLdouble)start[0];
  tmat[0][1] = (-sa); 
  tmat[1][1] = ca*cb; 
  tmat[2][1] = ca*sb; 
  tmat[3][1] = (GLdouble)start[1];
  tmat[0][2] = 0.0;   
  tmat[1][2] = (-sb); 
  tmat[2][2] = cb;    
  tmat[3][2] = (GLdouble)start[2];
  tmat[0][3] = 0.0;   
  tmat[1][3] = 0.0;   
  tmat[2][3] = 0.0;   
  tmat[3][3] = 1.0;

  glPushMatrix();
  glMultMatrixd((GLdouble*)tmat);
  gluCylinder( q, (GLfloat)radi, (GLfloat)0.0, height, n, 1);
  glPopMatrix();
}

void Model::ChangeArrows( const QString &new_arrow_type )
{
  arrow_type = new_arrow_type;
  panel->viewer->Redraw();
}

void Model::ChangeArrowColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  arrow_color[0] = color.red()/255.0;
  arrow_color[1] = color.green()/255.0;
  arrow_color[2] = color.blue()/255.0;
  panel->viewer->Redraw();
}

void Model::DrawArrows()
{
  glColor3f( arrow_color[0] , arrow_color[1] , arrow_color[2] );
  if ( ! (arrow_type == "None" ) )
  {
    if ( arrow_type == "Forces" )
    {

  double start[3],end[3];


  gluQuadricDrawStyle( q , GLU_FILL );
  glColor3d(0.00, 0.00, 1.00); // blue
  DrawCylinder( start , end , 0.1 , n );
}
    Atom *atom;
    ss = 10.0*panel->model->srs*scale;
    vb = vshift;
    vb *= ss;

    for ( iatom = 0 ; iatom < panel->geoeditor->mol->Natom(); iatom++ )
    {
      atom = panel->geoeditor->mol->atom[iatom];
      start[0] = atom->x;
      start[1] = atom->y;
      start[2] = atom->z;
      end[0] = atom->x + atom->force[0]*5.0;
      end[1] = atom->y + atom->force[1]*5.0;
      end[2] = atom->z + atom->force[2]*5.0;
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
    }
  };
}



