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
#ifndef GN_H
#define GN_H
 
#include "qapplication.h"
 
#include <Parameter.h>

#define SFB    0x00000010L           // Selecting for bond
#define SFA    0x00000020L           // Selecting for angle
#define SFD    0x00000040L           // Selecting for dihedral
#define SFRV   0x00000080L           // Selecting for distance value
#define SFAV   0x00000100L           // Selecting for angle value
#define SFDV   0x00000200L           // Selecting for diheral value
#define SFSURF 0x00000400L           // Selecting for surface
#define SFPG   0x00000800L           // Selecting for point group
#define SFFV   0x00001000L           // Selecting for fitting value
#define SFPB   0x00002000L           
#define QUEUED 0x00004000L           // GUI has been switched off     
 
class DwPanel;
class DwElement;
class System;

class GNagual : public QApplication 
{
  Q_OBJECT

  public:
    GNagual( System*, int argc = 0, char *argv[] = 0);

    char input_file[GN_MAX_STR_SIZE];
    long status;
 
    void SetCursor( const char* );
    void ErrorMessage(const char*, int, const char*, int);
    bool GetFileName(char*);
 
    System *sys;
    DwPanel *panel;

    DwElement *es;

  protected:

};
 
#endif // GN_H
