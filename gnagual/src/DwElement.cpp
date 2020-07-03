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
#include <QBoxLayout>
#include <QToolButton>

#include <Element.h>
#include <DwElement.h>

static int bpos[GN_MAX_ATOMIC_NUMBER+1][2] =  {{ 7,0},{ 0,0},{17,0},{ 0,1},
      { 1,1},{12,1},{13,1},{14,1},{15,1},{16,1},{17,1},{ 0,2},{ 1,2},
      {12,2},{13,2},{14,2},{15,2},{16,2},{17,2},{ 0,3},{ 1,3},{ 2,3},
      { 3,3},{ 4,3},{ 5,3},{ 6,3},{ 7,3},{ 8,3},{ 9,3},{10,3},{11,3},
      {12,3},{13,3},{14,3},{15,3},{16,3},{17,3},{ 0,4},{ 1,4},{ 2,4},
      { 3,4},{ 4,4},{ 5,4},{ 6,4},{ 7,4},{ 8,4},{ 9,4},{10,4},{11,4},
      {12,4},{13,4},{14,4},{15,4},{16,4},{17,4},{ 0,5},{ 1,5},{ 2,5},
      { 4,8},{ 5,8},{ 6,8},{ 7,8},{ 8,8},{ 9,8},{10,8},{11,8},{12,8},
      {13,8},{14,8},{15,8},{16,8},{17,8},{ 3,5},{ 4,5},{ 5,5},{ 6,5},
      { 7,5},{ 8,5},{ 9,5},{10,5},{11,5},{12,5},{13,5},{14,5},{15,5},
      {16,5},{17,5},{ 0,6},{ 1,6},{ 2,6},{ 4,9},{ 5,9},{ 6,9},{ 7,9},
      { 8,9},{ 9,9},{10,9},{11,9},{12,9},{13,9},{14,9},{15,9},{16,9},
      {17,9},{7,1},{7,2}}; 

DwElement::DwElement( QWidget *parent , int guess )
 //              : QDialog( parent ) 
{
  int i;

  selected_element = guess;
/*

  setWindowTitle( "Element Selector" );

  QGridLayout *mainLayout = new QGridLayout;

  for ( i = 0 ; i <= GN_MAX_ATOMIC_NUMBER ; i++ )
  {
    element[i] =  new QToolButton( this );
    element[i]->setText( ELEMENT_SYMBOL[i] ); 
    element[i]->setCheckable( true ); 
    connect(element[i],SIGNAL(clicked()),
            this,SLOT(SetDefaultElement())); 
  };
  element[selected_element]->setChecked( true );

  for ( i = 0 ; i <= GN_MAX_ATOMIC_NUMBER ; i++ )
  {
    mainLayout->addWidget( element[i] , bpos[i][1] , bpos[i][0] );
  };
  setLayout( mainLayout );
*/
}

void DwElement::SetDefaultElement() 
{
  int i; 

  for ( i = 0 ; i <= GN_MAX_ATOMIC_NUMBER ; i++ ) 
  {
    if ( element[i]->isChecked() && i != selected_element ) 
    {
      element[selected_element]->setChecked( false );
      selected_element = i;
    };
  };

  element[selected_element]->setChecked( true );

  accept();

}

int DwElement::GetSelectedElement( )
{
  return( selected_element );
}
