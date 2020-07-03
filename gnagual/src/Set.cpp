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
#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>

#include <Parameter.h>
#include <Set.h>
#include <Shell.h>
#include <Math.h>

using namespace std;

Set::Set(double ix,double iy,double iz)
{
  x = ix;
  y = iy;
  z = iz;
  shell.clear();
  nco = 0;
  is_neighbor = true;
};

void Set::LoadFromFile(char* basname,const char* atomname, bool potential)
{
  char line[GN_MAX_STR_SIZE];
  char filename[GN_MAX_STR_SIZE];

  const char* homeDir = getenv("HOME");

  if (potential)
  {
    if ((strncmp(basname,"ECP",3) == 0 )||
        (strncmp(basname,"QECP",4) == 0 )||
        (strncmp(basname,"RECP",4) == 0 )) 
      sprintf(filename,"%s/lib/potentials/ECP_PSEUDOS",homeDir);
    else
      sprintf(filename,"%s/lib/potentials/%s",homeDir,basname);
  }
  else
  {
    if ((strncmp(basname,"ECP",3) == 0 )||
        (strncmp(basname,"QECP",4) == 0 )||
        (strncmp(basname,"RECP",4) == 0 )) 
      sprintf(filename,"%s/lib/basis/ECP_BASIS",homeDir);
    else
      sprintf(filename,"%s/lib/basis/%s",homeDir,basname);
  }
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file [" << filename<< "]"<<endl;
    return;
  };
  string basisname = "(";
  basisname += basname;
  basisname += ")";
  string flag;
  if (potential) flag = "P-";
  else flag = "O-";
  flag += atomname;
  while(f.getline(line,GN_MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find(flag);
    if(sindex==0)
    {
      int tindex=stringLine.find(basisname);
      if (tindex>=0)
      {
        shell.clear();
        int k,n,iblk,nblk;
        double d,z;
        char jump[GN_MAX_STR_SIZE];
        f >> jump;
        if ( jump[0]=='#') 
        {
          f.getline(line,GN_MAX_STR_SIZE);
          f >> nblk;
        }
        else nblk = atoi(jump);
        int cp;
        for (iblk=0;iblk<nblk;iblk++)
        {
          Shell sb = Shell();
          f >> sb.n;
          f >> sb.l;
          f >> k;
          for (n=0;n<k;n++)
          {
            f >> z;
            f >> d;
            sb.z.push_back(z);
            sb.d.push_back(d);
          }
          if (!potential) 
          {
            sb.Normalize();
            sb.EvaluateRadius();
          }
          shell.push_back( new Shell( sb ) );
          cp = shell.size() - 1;
          if (cp==0) shell[cp]->ll = 0;
          else shell[cp]->ll = shell[cp-1]->ul + 1;
          shell[cp]->ul = shell[cp]->ll + shell[cp]->nco - 1;
        }
        nco = shell[cp]->ul + 1;
        break;
      }
    }
  }
  f.close();

  // Labeling 
  int i,lmax;
  if (potential) lmax = LMax();
  for (i=0;i<(signed)shell.size();i++)
  {
    if (potential)
    {
      if (shell[i]->l==lmax) shell[i]->type = CENTRAL_FORCE_SHELL;
      else shell[i]->type = SEMILOCAL_ECP_SHELL;
    }
    else
    {
      shell[i]->type = PRIMARY_BASIS_SHELL;
    }
  }

  return;
}

void Set::Print(char* filename)
{
  ofstream f(filename,ios::app);
  f <<endl;
  f << "Location ("<<x<<","<<y<<","<<z<<")"<<endl;
  f.close();
  for (int is=0;is<(signed)shell.size();is++)
  {
    shell[is]->Print(filename);
  }
}


int Set::LMax()
{
  int is,l,lmax;

  lmax = 0;
  for (is=0;is<(signed)shell.size();is++)
  {
    l = shell[is]->l;
    if (l>lmax) lmax = l;
  }
  return lmax;
}
