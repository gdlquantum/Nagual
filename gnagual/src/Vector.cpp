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
#include <stdio.h>

#include <iostream>

#include <Vector.h>
#include <Math.h>
#include <Parameter.h>

using namespace std;

// Constructors
Vector::Vector( double *v )
{
  x = v[0];
  y = v[1];
  z = v[2];
};

Vector::Vector( double xx , double yy , double zz ) 
{
  x = xx;
  y = yy;
  z = zz;
};

// Access
double& Vector::operator [] ( short n )
{
  if ( n == 0 )
  {
    return x;
  }
  else if ( n == 1 )
  {
    return y;
  }
  else if ( n == 2 )
  {
    return z;
  } 
  else
  {
    cout << "Asking for non existing element, returning z component"<<endl;
    return z;
  };
};

// this = result 

void Vector::operator += ( const Vector &cv )
{
  Vector v = cv;
  x += v[0];
  y += v[1];
  z += v[2];
};

void Vector::operator -= ( const Vector &cv )
{
  Vector v = cv; 
  x -= v[0];
  y -= v[1];
  z -= v[2];
};

void Vector::operator *= ( double s )
{
  x *= s;
  y *= s;
  z *= s;
};
void Vector::operator *= ( int n )
{
  x *= double(n);
  y *= double(n);
  z *= double(n);
};

void Vector::operator /= ( double s )
{
  x /= s;
  y /= s;
  z /= s;
};

void Vector::operator <= ( const Vector &cv )
{
  double m[3][3];
  Vector v = cv;

  m[0][0] = 0.0;   m[0][1] = -v[2]; m[0][2] = v[1];
  m[1][0] = v[2];  m[1][1] = 0.0;   m[1][2] = -v[0];
  m[2][0] = -v[1]; m[2][1] = v[0];  m[2][2] = 0.0;

  operator %= ( m );
};

void Vector::operator >= ( const Vector &v )
{
  operator <= ( v );
  operator *= ( -1.0 );
};

void Vector::operator %= ( double m[3][3] )
{
  double nx,ny,nz;

  nx = Dot( Vector( m[0][0] , m[0][1] , m[0][2] ) ); 
  ny = Dot( Vector( m[1][0] , m[1][1] , m[1][2] ) );
  nz = Dot( Vector( m[2][0] , m[2][1] , m[2][2] ) );

  x = nx;
  y = ny;
  z = nz;
};

void Vector::operator ^= ( double s ) 
{
  double norm = Norm();
  if ( norm > 0.0 )
  {
    operator *= ( s/norm );
  };
};


// this != result

Vector Vector::operator + ( const Vector &v )
{
  Vector u( x , y , z );
  u += v;
  return u;
};

Vector Vector::operator - ( const Vector &v )
{
  Vector u( x , y , z );
  u -= v;
  return u;
};

Vector Vector::operator < ( const Vector &v )
{
  Vector u( x , y , z );
  u <= v;
  return u;
};

Vector Vector::operator > ( const Vector &v )
{
  Vector u( x , y , z );
  u >= v;
  return u;
};

Vector Vector::operator % ( double m[3][3] )
{
  Vector u( x , y , z );
  u %= m;
  return u;
};

Vector Vector::operator ^ ( double s ) 
{
  Vector u( x , y , z );
  u ^= s;
  return u;
};

Vector Vector::operator * ( double s ) 
{
  Vector u( x , y , z );
  u *= s;
  return u;
};

Vector Vector::operator * ( int n ) 
{
  Vector u( x , y , z );
  u *= n;
  return u;
};

Vector Vector::operator / ( double s ) 
{
  Vector u( x , y , z );
  u /= s;
  return u;
};

// functions

double Vector::Dot( const Vector &cv )
{
  Vector v = cv;
  return( x*v[0] + y*v[1] + z*v[2] );
};

double Vector::Norm() 
{
  double r2 = x*x + y*y + z*z;
  if ( r2 < GN_TOL_NUM ) return 0.0;
  else return  sqrt(r2);
};

bool Vector::Colineal(const Vector &v)
{
  Vector u = v;

  if ( u.Norm() == 0.0 ) 
  {
    return true;
  }
  u ^= 1.0;
  u *= -Dot(u); 
  u.x += x;
  u.y += y;
  u.z += z;
  if (u.Norm()<GN_TOL_NUM) return true;
  else return false;
}

Vector Vector::Orthogonal() 
{
  double s;
  Vector u,v;

  if (Norm()<GN_TOL_NUM) return Vector( 1.0 , 0.0 , 0.0 );

  u = Vector( y , x , z );
  if (Colineal(u)) u = Vector( x , z , y );
  v = Vector( x , y , z );
  v ^= 1.0;
  s = u.Dot(v);
  v *= s;
  u -= v;
  u ^= 1.0;

  return u;
}

void Vector::Rotate( const Vector &v , double angle )
{
  double cosp,sinp,cost;
  double m[3][3];

  Vector u = v;
  if (Colineal(v)) return;

  if ( u.Norm() == 0.0 ) 
  {
    return;
  }
  u ^= 1.0;

  angle = angle*GN_PI/180.0;

  cosp = cos(angle);
  sinp = sin(angle); 
  cost = 1.0 - cosp;

  m[0][0] = u[0]*u[0]*cost + cosp;
  m[0][1] = u[0]*u[1]*cost - u[2]*sinp;
  m[0][2] = u[0]*u[2]*cost + u[1]*sinp;
  m[1][0] = u[1]*u[0]*cost + u[2]*sinp;
  m[1][1] = u[1]*u[1]*cost + cosp;
  m[1][2] = u[1]*u[2]*cost - u[0]*sinp;
  m[2][0] = u[2]*u[0]*cost - u[1]*sinp;
  m[2][1] = u[2]*u[1]*cost + u[0]*sinp;
  m[2][2] = u[2]*u[2]*cost + cosp;

  operator %= ( m );
};

// Transform Cartesian to Spherical coordinates. 
// History: - Creation (19.10.05, RFM)
//          - Translated to C++ (21.07.15, RFM)
void Vector::SphericalCoordinates(double *r, double *theta, double *phi)
{
  // Radius 
  *r = Norm();

  // Elevation angle 
  if (*r<GN_TOL_NUM) *theta = 0.0;
  else *theta = acos(z/(*r));

  // Azimuth angle 
  if (GN_ABS(x)<GN_TOL_NUM) 
  {
    if (GN_ABS(y)<GN_TOL_NUM) *phi = 0.0;
    else if (y<0.0) *phi = 1.5*GN_PI;
    else *phi = 0.5*GN_PI;
  }
  else if (x>0.0) *phi = atan(y/x);
  else *phi = atan(y/x) + GN_PI;
}

void Vector::Reflect( const Vector &normal )
{
  Vector other = Vector( normal );

  double nn = other.Dot(other);

  double rn = Dot(other);

  double factor = 2.0*rn/nn;
  other *= factor;
  x -= other.x;
  y -= other.y;
  z -= other.z;
  
}
