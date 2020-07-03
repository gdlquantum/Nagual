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
#ifndef GN_MATRIX_H
#define GN_MATRIX_H

#include <Matrix.h>

class Matrix
{
  public:
    Matrix(int,int);

    double *values;

    void operator += ( Matrix );
    double& operator () ( int, int );

    int NRow(void);
    int NCol(void);
    void Print(const char*,const char*,int,int);
    void Read(char*,char*,int,int);
    void Transpose(void);
    void Diagonalize(double*);
    void SVDPower(double*,double,double);
    void Symmetrize(void);
    void Fill(double);
    void SetZero(void);
    void SetValue( int, int , double );
    void ShiftValue( int, int , double );
    void GetRowValues( int, double* );
    void GetColValues( int, double* );
    void SetRowValues( int, double* );
    void SetColValues( int, double* );
    void Multiply(Matrix*,Matrix*,bool);
    void VectorMultiply(double*,double*);
    void Copy(Matrix*);
    void Add(Matrix*,int,int);
    void Scale(double);
    double Trace(void);
    double QTrace(Matrix*);

  protected:

    int dim_row;
    int dim_col;
};

#endif
