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
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <QtOpenGL/QGLWidget>
#include <QtCore>
#include <QtGui>
#include <QStringList>

#include <DwGeo.h>   
#include <GNagual.h>
#include <DwPanel.h>   
#include <Viewer.h>   
#include <Element.h>
#include <DwElement.h>
#include <Vector.h>
#include <Molecule.h>
#include <System.h>
#include <Atom.h>
#include <Units.h>

using namespace std;
using std::ifstream;


DwGeo::DwGeo( DwPanel *ip )
     : QWidget( 0 )
{
  panel = ip;
  gn = ip->gn;
  mol = panel->sys->mol;
  setWindowTitle( "Geometry Editor" );
  
  QList<QString> elebas;
  QList<QString> nucbas;

  QStringList Headers;
  Headers << "Atom" 
          << "Bond" << " " 
          << "Angle" << " " 
          << "Dihedral" << " " 
          << "e-Basis" << "Nuc-Basis";
  QGridLayout *mainLayout = new QGridLayout;

  tw = new QTableWidget( 1000 , 9 );
  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));
  tw->setHorizontalHeaderLabels(Headers);

  QToolButton *addButton =  new QToolButton( this );
  addButton->setText( "Add Atom" ); 
  connect(addButton,SIGNAL(clicked()),this,SLOT(AddAtomStart())); 

  QToolButton *delButton =  new QToolButton( this );
  delButton->setText( "Delete Last Atom" ); 
  connect(delButton,SIGNAL(clicked()),this,SLOT(DeleteAtom())); 

  QToolButton *fdButton =  new QToolButton( this );
  fdButton->setText( "Fit Distance" );
  connect(fdButton,SIGNAL(clicked()),this,SLOT(FitDistance()));

  timer = new QTimer(this);

  QComboBox *mComboBox = new QComboBox;
  mComboBox->addItem( QString("Measurement"));
  mComboBox->addItem( QString("Distance"));
  mComboBox->addItem( QString("Angle"));
  mComboBox->addItem( QString("Dihedral"));
  mComboBox->addItem( QString("Clean"));
  connect(mComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( SetupMeasurement(const QString &)));

  QCheckBox* cz = new QCheckBox( this );
  cz->setText( "Cartesians" ); 
  connect( cz , SIGNAL( clicked() ) , this , SLOT( ChangeCoordinates() ) );

  QCheckBox* cu = new QCheckBox( this );
  cu->setText( "Bohrs" ); 
  connect( cu , SIGNAL( clicked() ) , this , SLOT( ChangeUnits() ) );

  QComboBox *agComboBox = new QComboBox;
  agComboBox->addItem( QString("Add Group"));
  agComboBox->addItem( QString("--OH"));
  agComboBox->addItem( QString("--NH2"));
  agComboBox->addItem( QString("--CH3"));
  agComboBox->addItem( QString("File XYZ"));
  connect(agComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( AddGroupStart(const QString &)));


  QComboBox *pgComboBox = new QComboBox;
  pgComboBox->addItem( QString("Symmetry Group"));
  pgComboBox->addItem( QString("Td"));
  pgComboBox->addItem( QString("Oh"));
  pgComboBox->addItem( QString("Ih"));
  pgComboBox->addItem( QString("Cs"));
  pgComboBox->addItem( QString("Ci"));
  pgComboBox->addItem( QString("C1"));
  pgComboBox->addItem( QString("C2h"));
  pgComboBox->addItem( QString("C3h"));
  pgComboBox->addItem( QString("C4h"));
  pgComboBox->addItem( QString("C5h"));
  pgComboBox->addItem( QString("C2v"));
  pgComboBox->addItem( QString("C3v"));
  pgComboBox->addItem( QString("C4v"));
  pgComboBox->addItem( QString("C5v"));
  pgComboBox->addItem( QString("C6v"));
  pgComboBox->addItem( QString("C7v"));
  pgComboBox->addItem( QString("C8v"));
  pgComboBox->addItem( QString("C9v"));
  pgComboBox->addItem( QString("C*v"));
  pgComboBox->addItem( QString("S2"));
  pgComboBox->addItem( QString("S4"));
  pgComboBox->addItem( QString("S6"));
  pgComboBox->addItem( QString("S8"));
  pgComboBox->addItem( QString("D2h"));
  pgComboBox->addItem( QString("D3h"));
  pgComboBox->addItem( QString("D4h"));
  pgComboBox->addItem( QString("D5h"));
  pgComboBox->addItem( QString("D6h"));
  pgComboBox->addItem( QString("D7h"));
  pgComboBox->addItem( QString("D8h"));
  pgComboBox->addItem( QString("D9h"));
  pgComboBox->addItem( QString("D*h"));
  pgComboBox->addItem( QString("D3d"));
  pgComboBox->addItem( QString("D5d"));
  pgComboBox->addItem( QString("D7d"));
  pgComboBox->addItem( QString("D9d"));
  pgComboBox->addItem( QString("s1c1"));
  connect(pgComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( BuildGroup(const QString &)));

  mainLayout->addWidget( addButton   , 0 , 0 , 1, 1 );
  mainLayout->addWidget( delButton   , 0 , 1 , 1, 1 );
  mainLayout->addWidget( agComboBox  , 1 , 0 , 1, 1 );
  mainLayout->addWidget( pgComboBox  , 1 , 1 , 1, 1 );
  mainLayout->addWidget( cz          , 1 , 2 , 1, 1 );
  mainLayout->addWidget( fdButton    , 2 , 0 , 1, 1 );
  mainLayout->addWidget( mComboBox   , 2 , 1 , 1, 1 );
  mainLayout->addWidget( cu          , 2 , 2 , 1, 1 );
  mainLayout->addWidget( tw          , 4 , 0 , 3, -1 );
  setLayout( mainLayout );

  CCSET = ZMatrix;
  default_element = 6;
  use_default = false;
  bohr = false;
  strcpy(group,"Atom");
}

void DwGeo::AddAtomStart( )
{
  if ( !use_default )
  {
    //gn->es->exec(); 
    default_element = gn->es->GetSelectedElement();
  };
  Enabled( false );
  if ( mol->Natom() < 2 )
  {
    AddAtom( mol->Natom(),0, 0 );
  }
  else 
  {
    if ( mol->Natom() < 3 )
    {
      gn->panel->Report( "Select one atom" ); 
    }
    else if ( mol->Natom() < 4 )
    {
      gn->panel->Report( "Select two atoms" );
    }
    else
    {
      gn->panel->Report( "Select three atoms" );  
    };
    Enabled( false );
    gn->status |= SFB;
    gn->SetCursor( "Selection" );
  };
}

void DwGeo::Update()
{
  int iatom;
  Atom *a;

  disconnect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
             this , SLOT( Edited( QTableWidgetItem*) ));

  for ( iatom = 0 ; iatom < mol->Natom() ; iatom++ )
  {
    a = mol->GetAtom(iatom);
    if ( ! tw->item(iatom,0) )
    {
      tw->setItem(iatom,0, new QTableWidgetItem( " " ));
      tw->setItem(iatom,1, new QTableWidgetItem( " " ));
      tw->setItem(iatom,2, new QTableWidgetItem( " " ));
      tw->setItem(iatom,3, new QTableWidgetItem( " " ));
      tw->setItem(iatom,4, new QTableWidgetItem( " " ));
      tw->setItem(iatom,5, new QTableWidgetItem( " " ));
      tw->setItem(iatom,6, new QTableWidgetItem( " " ));
    }
    (tw->item(iatom,0))->setText(QString( ELEMENT_SYMBOL[a->atomic_number] ) );
    if (CCSET==ZMatrix) 
    {
      if ( iatom > 0 ) 
      {
        // Remember: users sees angstroms unless he asks otherwise, 
        // inside we have bohrs always
        (tw->item(iatom,1))->setText( QString::number( a->ref_bond ) );
        if (bohr) (tw->item(iatom,2))->setText( QString::number( a->bond ) );
        else (tw->item(iatom,2))->setText( QString::number( BohrToAngstrom(a->bond) ) );
        if ( iatom > 1 )
        {
          (tw->item(iatom,3))->setText( QString::number( a->ref_angle ) );
          (tw->item(iatom,4))->setText( QString::number( a->angle ) );
          if ( iatom > 2 )
          {
            (tw->item(iatom,5))->setText(QString::number( a->ref_dihedral ) );
            (tw->item(iatom,6))->setText(QString::number( a-> dihedral ) );
          }
          else
          {
            (tw->item( iatom , 5 ))->setText( " " ); 
            (tw->item( iatom , 6 ))->setText( " " );
          }
        }
        else
        {
          (tw->item( iatom , 3 ))->setText( " " ); 
          (tw->item( iatom , 4 ))->setText( " " );
          (tw->item( iatom , 5 ))->setText( " " );
          (tw->item( iatom , 6 ))->setText( " " );
        }
      }
      else
      {
        (tw->item( iatom , 1 ))->setText( " " ); 
        (tw->item( iatom , 2 ))->setText( " " );
        (tw->item( iatom , 3 ))->setText( " " );
        (tw->item( iatom , 4 ))->setText( " " );
        (tw->item( iatom , 5 ))->setText( " " );
        (tw->item( iatom , 6 ))->setText( " " );
      }
    }
    else
    { 
      (tw->item( iatom , 1 ))->setText( QString( "x" ) );
      (tw->item( iatom , 3 ))->setText( QString( "y" ) );
      (tw->item( iatom , 5 ))->setText( QString( "z" ) );
      // Remember: users sees angstroms unless he asks otherwise, 
      // inside we have bohrs always
      if (bohr) 
      {
        (tw->item( iatom , 2 ))->setText( QString::number( a->x ) );
        (tw->item( iatom , 4 ))->setText( QString::number( a->y ) );
        (tw->item( iatom , 6 ))->setText( QString::number( a->z ) );
      }
      else
      {
        (tw->item(iatom,2))->setText( QString::number(BohrToAngstrom(a->x)) );
        (tw->item(iatom,4))->setText( QString::number(BohrToAngstrom(a->y)) );
        (tw->item(iatom,6))->setText( QString::number(BohrToAngstrom(a->z)) );
      }
    }
  }

  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));

}

void DwGeo::DeleteAtom()
{
  int natom = mol->Natom();
  if ( natom > 0 )
  {
    mol->DeleteLastAtom();
    natom--;
    (tw->item( natom , 0 ))->setText( QString( " " ) );
    (tw->item( natom , 1 ))->setText( QString( " " ) );
    (tw->item( natom , 2 ))->setText( QString( " " ) );
    (tw->item( natom , 3 ))->setText( QString( " " ) );
    (tw->item( natom , 4 ))->setText( QString( " " ) );
    (tw->item( natom , 5 ))->setText( QString( " " ) );
    (tw->item( natom , 6 ))->setText( QString( " " ) );
    panel->viewer->Redraw();
    Enabled( true );
  };
}

void DwGeo::Enabled( bool enabled )
{
  setEnabled( enabled );
}

void DwGeo::Edited( QTableWidgetItem* item ) 
{
  int iatom = tw->currentRow();
  int col = tw->currentColumn();

  if ( col == 0 )
  {
    mol->atom[iatom]->SetSymbol( (item->text().toLatin1()).data() );
    item->setText( QString( ELEMENT_SYMBOL[mol->atom[iatom]->atomic_number] ) );
  }
  else if ( col == 1 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_bond ) );
  }
  else if ( col == 2 )
  {
    if (CCSET==ZMatrix) 
    {
      if (bohr) mol->atom[iatom]->bond = item->text().toDouble();
      else mol->atom[iatom]->bond = AngstromToBohr( item->text().toDouble() );
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->x = item->text().toDouble();
      else mol->atom[iatom]->x = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  }
  else if ( col == 3 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_angle ) );
  }
  else if ( col == 4 )
  {
    if (CCSET==ZMatrix) 
    {
      mol->atom[iatom]->angle = item->text().toDouble();
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->y = item->text().toDouble();
      else mol->atom[iatom]->y = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  }
  else if ( col == 5 )
  {
    item->setText( QString::number( mol->atom[iatom]->ref_dihedral ) );
  }
  else if ( col == 6 )
  {
    if (CCSET==ZMatrix) 
    {
      mol->atom[iatom]->dihedral = item->text().toDouble();
      mol->Z2C(iatom);
    }
    else 
    {
      if (bohr) mol->atom[iatom]->z = item->text().toDouble();
      else mol->atom[iatom]->z = AngstromToBohr( item->text().toDouble() );
      mol->C2Z(iatom);
    }
  } 
  panel->viewer->Redraw();
}

void DwGeo::BuildGroup( const QString &s )
{
  if ( s == "Symmetry Group" ) return;
  sprintf( gname , "%s", (s.toLatin1()).data() );
  gn->status |= SFPG;
  Enabled ( false );
}

void DwGeo::AddGroupStart(const QString &s)
{
  bool hold;

  strcpy(group, (s.toLatin1()).data() );

  hold = use_default;
  use_default = true;

  AddAtomStart();

  use_default = hold;
}


void DwGeo::ChangeCoordinates()
{
  if (CCSET==Cartesians)
  {
    CCSET = ZMatrix;
  }
  else
  {
    CCSET = Cartesians;
  }
  Update();
}

void DwGeo::UseDefault()
{
  use_default =  (use_default ? false : true );
}

void DwGeo::FitDistance( void )
{
  gn->SetCursor( "Selection" );
  gn->status |= SFFV;
  gn->panel->Report("Select two atoms");
}

int DwGeo::SGDrv(int* ref) 
{
  if ( gname[0] == 's' )
  {
    return mol->SGBDrv(gname,ref);
  }
  else 
  {
    return mol->PGBDrv(gname,ref);
  };
};

void DwGeo::ChangeUnits()
{
  if (bohr)
  {
    bohr = false;
  }
  else
  {
    bohr = true;
  }
  Update();
}

// Reading from files is mediated by windows
// (File substitutes user, not window)
void DwGeo::Read(char* filename, char* fmt)
{
  if ( strcasecmp( fmt , "XYZ" ) == 0 ) mol->ReadXYZ(filename);
  else if ( strcasecmp( fmt , "ZMT" ) == 0 ) mol->ReadZMatrix(filename);
  Update(); 
  panel->viewer->Redraw();
}


void DwGeo::SetupMeasurement( const QString &s )
{
  if ( s == "Measurement" ) return;
  if ( s == "Distance" && mol->Natom() > 1 ) 
  {
    gn->SetCursor( "Selection" );
    gn->status |= SFRV;
    panel->Report("Select two atoms");
  }
  else if ( s == "Angle" && mol->Natom() > 2)
  {
    gn->SetCursor( "Selection" );
    gn->status |= SFAV;
    panel->Report("Select three atoms");
  }
  else if ( s == "Dihedral" && mol->Natom() > 3 )
  {
    gn->SetCursor( "Selection" );
    gn->status |= SFDV;
    panel->Report("Select four atoms");
  }
  else if ( s == "Clean" && panel->viewer ) 
  {
    panel->viewer->monitored_bonds.clear();
    panel->viewer->monitored_angles.clear();
    panel->viewer->monitored_dihedrals.clear();
    panel->viewer->Redraw();
  };
}

void DwGeo::AddAtom( int na, int nb, int nc )
{
  double bond;

  if ( strncmp( group , "Atom", 4) == 0 )
  {
    if (na>0) bond = ELEMENT_COV_R[default_element] + 
                   ELEMENT_COV_R[mol->atom[na-1]->atomic_number];
    else bond = 0.0;
    mol->AddAtom( default_element, na - 1, nb - 1, nc - 1,
                         bond,109.467, 180.0);
  }
  else
  {
    AddGroup(na-1,nb-1,nc-1);
    strcpy( group , "Atom");
  }
  Update();
  panel->viewer->Redraw();
  Enabled( true );
}

void DwGeo::AddGroup( int n1, int n2, int n3 )
{
  int an,na,nb,nc;
  double bond;

  if ( strncmp(group,"--OH",4) == 0 ||
       strncmp(group,"--NH",4) == 0 ||
       strncmp(group,"--CH",4) == 0 ) 
  {   
    if ( strncmp(group,"--OH",4) == 0 ) an = 8;
    else if ( strncmp(group,"--NH",4) == 0 ) an = 7;
    else if ( strncmp(group,"--CH",4) == 0 ) an = 6;

    na = n1;
    nb = n2;
    nc = n3;
    mol->atom[na]->SetAtomicNumber( an );
    bond = ELEMENT_COV_R[an] + ELEMENT_COV_R[1];
    mol->AddAtom(1,na,nb,nc,bond,109.467,180.0);
    if ( strncmp(group,"--OH",4) != 0 ) 
    {
      mol->AddAtom(1,na,nb,nc,bond,109.467,300.0);
      if ( strncmp(group,"--CH",4) == 0 ) 
        mol->AddAtom(1,na,nb,nc,bond,109.467,60.0);
    }
  }
  else if ( strncmp(group,"File XYZ",8) == 0 || 
            strncmp(group,"File ZMT",8) == 0 )
  {
    char filename[GN_MAX_STR_SIZE];
    int iatom,patom;
    double angle,dihedral;

    patom = mol->Natom() - 1;

    Molecule mol2 = Molecule();

    if (!gn->GetFileName(filename))
      cout << "Unable to open file "<<endl;

    if ( strncmp(group,"File XYZ",8) == 0 ) mol2.ReadXYZ(filename);
    else if ( strncmp(group,"File ZMT",8) == 0 ) mol2.ReadZMatrix(filename);

    for (iatom=0;iatom<mol2.Natom();iatom++)
    {
      an = mol2.atom[iatom]->atomic_number;
      if (iatom==0)  
      {
        mol->atom[n1]->SetAtomicNumber( an );
      }
      else
      {
        na = patom + mol2.atom[iatom]->ref_bond;
        if (na==patom) na = n1;
        bond = mol2.atom[iatom]->bond;

        if (iatom<2)
        {
          nb = n2;
          angle = 109.467;
        }
        else
        {
          nb = patom + mol2.atom[iatom]->ref_angle;
          angle = mol2.atom[iatom]->angle;
        }
        if (nb==patom) nb = n1;

        if (iatom<3)
        {
          nc = n3;
          dihedral = 60.0 + 120.0*iatom;
        }
        else
        {
          nc = patom + mol2.atom[iatom]->ref_dihedral;
          dihedral = mol2.atom[iatom]->dihedral;
        }
        if (nc==patom) nc = n1;

        mol->AddAtom(an,na,nb,nc,bond,angle,dihedral);
      }
    }
  }
}
