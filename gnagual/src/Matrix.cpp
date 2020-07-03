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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Matrix.h>
#include <Math.h>

using namespace std;

Matrix::Matrix(int n,int m)
{
  dim_row = n;
  dim_col = m;

  values = new double[dim_row*dim_col];
  SetZero();
}

double & Matrix::operator () ( int i, int j )
{
  return values[i*dim_col+j];
}

void Matrix::operator += ( Matrix other )
{
  int i,j,k;

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    for (j=0;j<dim_col;j++)
    {
      values[k] += other(i,j);
      k++;
    }
  }
}


void Matrix::Print(const char* filename, const char* header, int w, int nd)
{

  ofstream f(filename,ios::app);
  f << "=== "<<header<<" ==="<<endl;

  int nblk = 75/w;
  int ll,ul;
  ll = 0;
  while (ll<dim_col)
  {
    ul = GN_MIN(dim_col-1,ll+nblk-1);
    f << "     ";
    for (int j=ll;j<=ul;j++)
      f << setw(w) << setprecision(nd) << fixed << j;
    f << endl;
    for (int i=0;i<dim_row;i++)
    {
      f << setw(4) << i << " ";
      for (int j=ll;j<=ul;j++)
        f << setw(w) << setprecision(nd) << fixed << values[i*dim_col+j];
      f << endl;
    }
    f << endl;
    ll = ul + 1;
  }
  f.close();
}

int Matrix::NRow()
{
  return dim_row;
}
int Matrix::NCol()
{
  return dim_col;
}

void Matrix::GetRowValues( int id, double *val )
{
  int ll = id*dim_col;
  for (int j=0;j<dim_col;j++)
    val[j] = values[ll+j];
}
void Matrix::GetColValues( int id, double *val )
{
  for (int j=0;j<dim_row;j++)
    val[j] = values[j*dim_col+id];
}
void Matrix::SetRowValues( int id, double *val )
{
  int ll = id*dim_col;
  for (int j=0;j<dim_col;j++)
    values[ll+j] = val[j];
}
void Matrix::SetColValues( int id, double *val )
{
  for (int j=0;j<dim_row;j++)
    values[j*dim_col+id] = val[j];
}

void Matrix::SetValue( int i, int j, double val)
{
  values[i*dim_col+j] = val;
}

void Matrix::ShiftValue( int i, int j, double val)
{
  values[i*dim_col+j] += val;
}

// Notice that transposition switches dimensionality
// Roberto Flores-Moreno, Ago 2015
void Matrix::Transpose(void)
{
  int i,j;
  int dim_row_h,dim_col_h;
  double *h;

  dim_row_h = dim_col;
  dim_col_h = dim_row;
  h = new double[dim_col_h*dim_row_h];

  for (i=0;i<dim_row;i++)
    for (j=0;j<dim_col;j++)
      h[j*dim_col_h+i] = values[i*dim_col+j];

  dim_row = dim_row_h;
  dim_col = dim_col_h;
  for (i=0;i<dim_row;i++)
    for (j=0;j<dim_col;j++)
      values[i*dim_col+j] = h[i*dim_col+j];

  delete[] h;
}

void Matrix::Multiply(Matrix *other,Matrix *res,bool right)
{
  int i,j,k;
  double s;

  if (right)
  { 
    if (dim_col != other->NRow() ) 
    {
      cout << "WARNING: Matrix dimensions are diferent" <<endl;
    }
    double others[dim_col];
    for (i=0;i<dim_row;i++)
    {
      for (j=0;j<other->NCol();j++)
      {
        other->GetColValues(j,others);
        s = 0.0;
        for (k=0;k<dim_col;k++)
        {
          s += values[i*dim_col+k]*others[k];
        } 
        res->SetValue(i,j,s);
      }
    }
  }
  else
  {
    if (other->NCol() != dim_row ) 
    {
      cout << "WARNING: Matrix dimensions are diferent" <<endl;
    }
    double others[dim_row];
    for (i=0;i<other->NRow();i++)
    {
      other->GetRowValues(i,others);
      for (j=0;j<dim_col;j++)
      {
        s = 0.0;
        for (k=0;k<dim_row;k++)
        {
          s += others[k]*values[k*dim_col+j];
        } 
        res->SetValue(i,j,s);
      }
    }
  }
}

void Matrix::Fill(double val)
{
  int i,j,k;

  k = 0;
  for (i=0;i<dim_row;i++)
  {
    for (j=0;j<dim_col;j++)
    {
      values[k] = val;
      k++;
    }
  }
}

void Matrix::SetZero() { Fill(0.0); }

void Matrix::Symmetrize(void)
{
  int i,j;

  for (i=0;i<dim_row;i++)
    for (j=i+1;j<dim_col;j++)
      values[j*dim_col+i] = values[i*dim_col+j];
}

double Matrix::Trace()
{
  int i;
  double t;

  t = 0.0;
  for (i=0;i<dim_row;i++)
    t += values[i*dim_col+i];
  return t;
}

double Matrix::QTrace(Matrix *other)
{
  int i,j,ll;
  double t;
  double others[dim_col];

  if ((other->NRow()!=dim_row)||(other->NCol()!=dim_col))
  {
    cout << "Attemp to trace matrices of different dimensions "<<endl;
    exit(EXIT_FAILURE);
  }

  t = 0.0;
  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i,others);
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      t += others[j]*values[ll+j];
    }
  }
  return t;
}

void Matrix::Copy(Matrix *other)
{
  int i,j,ll;
  double others[dim_col];

  for (i=0;i<dim_row;i++)
  {
    other->GetRowValues(i,others);
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      values[ll+j] = others[j];
    }
  }
}

void Matrix::VectorMultiply(double *in,double *out)
{
  int i,ii,k;
  double s;

  for (i=0;i<dim_row;i++)
  {
    ii = i*dim_col;
    s = 0.0;
    for (k=0;k<dim_col;k++)
      s += values[ii+k]*in[k];
    out[i] = s;
  }
}

void Matrix::Add(Matrix *other, int llrow, int llcol)
{
  int i,j,ll;
  double others[other->NCol()];

  for (i=llrow;i<llrow+other->NRow();i++)
  {
    other->GetRowValues(i-llrow,others);
    ll = i*dim_col;
    for (j=llcol;j<llcol+other->NCol();j++)
    {
      values[ll+j] += others[j-llcol];
    }
  }
}

void Matrix::Scale(double s)
{
  int i,j,ll;

  for (i=0;i<dim_row;i++)
  {
    ll = i*dim_col;
    for (j=0;j<dim_col;j++)
    {
      values[ll+j] *= s;
    }
  }
}

// Loading matrix from file (AVR, Jul 2015)
// Generalized to any matrix by RFM, Aug 2015
// Adequate for print width of 20
void Matrix::Read(char *filename,char *header, int id, int ncol)
{
  int jd, i, j, k, ll,ul,matfrac, resid, a; //

  SetZero();

  resid = (ncol)%3; // Tells how many rows are in the last matrix
  matfrac = (ncol)/3+1; // Tells in how many blocks the matrix was divided
  string line;

  ifstream file (filename, ios:: in);

  if(file.is_open())
  {
    jd = -1;
    while (jd!=id)
    {
      getline(file,line);
      if ( strncmp(line.c_str(),"Subsystem ID:",13) == 0 )
        jd = atoi(&line[13]);
    }
    for ( ; line != header;) // Ignore lines
      getline(file,line);
    for (k = 0; k < matfrac-1; k++) // Counts matrix parts
    {
      ll = 3*k;
      if (k==matfrac-1) ul = 3*k + resid;
      else ul = 3*(k+1);
      file >> a; file >> a; file >> a;
      for (i = 0; i < dim_row; i++)
      {
        file >> a;
        for(j = ll; j < ul; j++)  // Counts 3 columns and ends loop
          file >> values[i*dim_col+j];
      }
    }
  }
  else cout << "Unable to open file " << filename << endl;
  
  file.close();
}

