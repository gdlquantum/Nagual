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
#ifndef GN_VECTOR_H
#define GN_VECTOR_H

class Vector
{
public:
  Vector( double xx = 0.0 ,
          double yy = 0.0 ,
          double zz = 0.0 );
  Vector( double* );

  double x;
  double y;
  double z; 

  Vector operator + ( const Vector & );
  Vector operator - ( const Vector & );
  Vector operator * ( double );
  Vector operator * ( int );
  Vector operator / ( double );
  Vector operator % ( double m[3][3] );
  Vector operator ^ ( double );
  Vector operator < ( const Vector & );
  Vector operator > ( const Vector & );

  void operator += ( const Vector & );
  void operator -= ( const Vector & );
  void operator >= ( const Vector & );
  void operator <= ( const Vector & );
  void operator *= ( double );
  void operator *= ( int );
  void operator /= ( double );
  void operator %= ( double m[3][3] );
  void operator ^= ( double );

  double& operator [] ( short n );

  bool Colineal(const Vector &);
  double Dot( const Vector & );
  double Norm(void);
  void SphericalCoordinates(double*,double*,double*);
  void Rotate( const Vector & , double );
  void Rotate2( const Vector & , double );
  void Reflect( const Vector & );
  Vector Orthogonal(void);

};

#endif

