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
#include <string.h>

#include <iostream>
#include <fstream>

#include <System.h>
#include <Spin.h>
#include <Set.h>
#include <Shell.h>
#include <Molecule.h>
#include <Element.h>
#include <Atom.h>
#include <Units.h>
#include <Math.h>
#include <Matrix.h>

using namespace std;

System::System()
{
  mol = new Molecule();
  nspin = 0;
  charges = 0;

  temperature = 0.0;    // Standard orbital calculations
  chemical_potential = 0.0;   
};

void System::Setup(QChem *iqchem)
{
  qchem = iqchem;
}

void System::Print(char* filename)
{

};

double System::GetAtomCharge(int id)
{
  double ret;

  ret = 0.0;

  if (id>=0&&id<mol->Natom())
    if (charges) ret = charges[id];

  return ret;
}

double System::Number(int type, int spin) 
{ 
// return number[type][spin]; 
} 

void System::ChangeNumber( double n, int type, int spin ) 
{ 
// number[type][spin] = n; 
} 

void System::AddSpin(char* name,double mass,double charge,int number_of_particles)
{
  spin[nspin] = new Spin(this,nspin,name,mass,charge,number_of_particles);
  nspin++;
}

void System::Read(char* filename,char* fmt)
{
  if ( strcasecmp( fmt , "NGL" ) == 0 ) ReadNGL(filename);
  else if ( strcasecmp( fmt , "XYZ" ) == 0 ) mol->ReadXYZ(filename);
}

void System::ReadNGL(char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    cout << "Error opening file "<<filename<<endl;
    return;
  };

  mol->Clear();

  int natom,nspin;
  char line[GN_MAX_STR_SIZE];
  char str[GN_MAX_STR_SIZE];
  while (f.getline(line,GN_MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find("Number of atoms:"); 
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    { 
      int i,iatom,iblk,ispin,k,l,n,nblk;
      double d,x,y,z;
      Atom *a;

      natom = atoi(&line[17]);
      f.getline(line,GN_MAX_STR_SIZE);
      nspin = atoi(&line[25]);
      f.getline(line,GN_MAX_STR_SIZE);
      f.getline(line,GN_MAX_STR_SIZE);
      for (i=0;i<natom;i++) 
      {
        f >> str;
        f >> x;
        f >> y;
        f >> z;
        x = AngstromToBohr(x);
        y = AngstromToBohr(y);
        z = AngstromToBohr(z);
        a = new Atom(str,x,y,z);
        mol->atom.push_back( a );
      }
      mol->Center();
      f.getline(line,GN_MAX_STR_SIZE);
      f.getline(line,GN_MAX_STR_SIZE);
      f.getline(line,GN_MAX_STR_SIZE);
      for (ispin=0;ispin<nspin;ispin++) 
      {
        char name[GN_MAX_STR_SIZE];
        int npart;
        double mass,charge;

        f >> str;
        f >> str; 
        f >> name; 
        f >> str;
        f >> npart; 
        f >> str;
        f >> mass; 
        f >> str;
        f >> charge; 
        AddSpin(name,mass,charge,npart);
        spin[ispin]->basis.clear();
        spin[ispin]->nco = 0;
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        for (iatom=0;iatom<natom;iatom++) 
        {
          Set *set;

          set = new Set( mol->atom[iatom]->x, 
                          mol->atom[iatom]->y, 
                          mol->atom[iatom]->z );
          if (iatom==0) set->ll = 0;
          else set->ll = spin[ispin]->basis[iatom-1]->ul + 1;
          set->nco = 0;
          f.getline(line,GN_MAX_STR_SIZE);
          f.getline(line,GN_MAX_STR_SIZE);  // %%%
          f.getline(line,GN_MAX_STR_SIZE);
          f.getline(line,GN_MAX_STR_SIZE);
          f >> nblk;
          for (iblk=0;iblk<nblk;iblk++) 
          {
            Shell *shell;

            shell = new Shell();
            if (iblk==0) shell->ll = 0;
            else shell->ll = set->shell[iblk-1]->ul + 1;
            f >> shell->n;
            f >> shell->l;
            shell->nco = ((shell->l+1)*(shell->l+2))/2;
            shell->ul = shell->ll + shell->nco - 1;
            set->nco = set->nco + shell->nco;
            f >> k;
            shell->z.clear();
            shell->d.clear();
            for (i=0;i<k;i++) 
            {
              f >> z;
              f >> d;
              shell->z.push_back(z);
              shell->d.push_back(d);
            }
            shell->Normalize();
            set->shell.push_back( shell );
          }
          set->ul = set->ll + set->nco - 1;
          spin[ispin]->basis.push_back( set );
          spin[ispin]->nco = spin[ispin]->nco + set->nco;
          f.getline(line,GN_MAX_STR_SIZE);
        }
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
      }
    } 
    int charindex=stringLine.find("Population analysis"); 
    if((charindex>=0)&&(charindex<(signed)stringLine.size())) 
    { 
      f.getline(line,GN_MAX_STR_SIZE);
      f.getline(line,GN_MAX_STR_SIZE);
      for (int i=0;i<mol->Natom();i++) 
      {
        f >> str;
        f >> mol->atom[i]->charge;
      }
    }
    int oindex=stringLine.find("Orbital energies and coefficients"); 
    if((oindex>=0)&&(oindex<(signed)stringLine.size())) 
    { 
      int i,ibatch,ispin,j,nbatch,nco,norb,top;
      double val[5];

      for (ispin=0;ispin<nspin;ispin++) 
      {
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f >> str;
        f >> str;
        f >> str;
        f >> norb;
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        spin[ispin]->oo = new double[norb];
        for (i=0;i<norb;i++)
        {
          f >> str;
          f >> spin[ispin]->oo[i];
          f.getline(line,GN_MAX_STR_SIZE);
        }
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        f.getline(line,GN_MAX_STR_SIZE);
        nbatch = norb/5;
        if (GN_MOD(norb,5)!=0) nbatch = nbatch + 1;
        nco = spin[ispin]->nco;
        spin[ispin]->C = new Matrix(nco,nco);
        for (ibatch=0;ibatch<nbatch;ibatch++)
        {
          if (ibatch==nbatch-1) 
          {
            if (GN_MOD(norb,5)==0) top = 5;
            else top = norb - 5*(norb/5);
          }
          else top = 5;
          for (j=0;j<top;j++)
            f >> val[j];
          for (i=0;i<nco;i++)
          {
            for (j=0;j<top;j++)
              f >> val[j];
            for (j=0;j<top;j++)
              spin[ispin]->C->SetValue(ibatch*5+j,i,val[j]);
          }
        }
        spin[ispin]->P.push_back( new Matrix(nco,nco) );
        spin[ispin]->BuildDensityMatrix();
      }
    }
    int gindex=stringLine.find("Molecular forces"); 
    if((gindex>=0)&&(gindex<(signed)stringLine.size())) 
    {
      f.getline(line,GN_MAX_STR_SIZE);
      for (int iatom=0;iatom<natom;iatom++) 
      {
        f >> str;
        f >> mol->atom[iatom]->force[0];
        f >> mol->atom[iatom]->force[1];
        f >> mol->atom[iatom]->force[2];
      }
    }
    int hindex=stringLine.find("Responses to density matrices"); 
    if((hindex>=0)&&(hindex<(signed)stringLine.size())) 
    { 
      int ispin,ipert,nco,npert;

      f.getline(line,GN_MAX_STR_SIZE);
      f >> str;
      f >> str;
      f >> str;
      f >> npert;
      for (ispin=0;ispin<nspin;ispin++) 
      {
        nco = spin[ispin]->nco;
        for (ipert=0;ipert<npert;ipert++) 
          spin[ispin]->P.push_back( new Matrix(nco,nco) );
      }
    }
    int pindex=stringLine.find("Density matrix response"); 
    if((pindex>=0)&&(pindex<(signed)stringLine.size())) 
    { 
      int i,ibatch,ipert,ispin,j,nbatch,nco,norb,top;
      double val[5];

      f.getline(line,GN_MAX_STR_SIZE);
      f >> str;
      f >> ispin;
      ispin = ispin - 1;
      f >> str;
      f >> ipert;
      f.getline(line,GN_MAX_STR_SIZE);
      f.getline(line,GN_MAX_STR_SIZE);
      nco = spin[ispin]->nco;
      nbatch = nco/5;
      if (GN_MOD(nco,5)!=0) nbatch = nbatch + 1;
      for (ibatch=0;ibatch<nbatch;ibatch++)
      {
        if (ibatch==nbatch-1) 
        {
          if (GN_MOD(nco,5)==0) top = 5;
          else top = nco - 5*(nco/5);
        }
        else top = 5;
        for (j=0;j<top;j++)
          f >> str;
        for (i=0;i<nco;i++)
        {
          for (j=0;j<top;j++)
            f >> val[j];
          for (j=0;j<top;j++)
            spin[ispin]->P[ipert]->SetValue(i,ibatch*5+j,val[j]);
        }
      }
    }
  }
  f.close();
}


