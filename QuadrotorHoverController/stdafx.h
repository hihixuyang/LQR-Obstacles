// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define _CRT_RAND_S
#define _USE_MATH_DEFINES

#include "targetver.h"

#include <math.h>
#include <stdio.h>
#include <tchar.h>
#include <vector>
#include <deque>
#include "matrix.h"
#include "callisto.h"
#include "callistoTypes.h"
#include "time.h"


inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) {
  Matrix<4,1> q;

  q[0] = 0.5*sqrt(std::max(1+R(0,0)-R(1,1)-R(2,2),0.0))*(R(2,1) - R(1,2) >= 0 ? 1 : -1);
  q[1] = 0.5*sqrt(std::max(1-R(0,0)+R(1,1)-R(2,2),0.0))*(R(0,2) - R(2,0) >= 0 ? 1 : -1);
  q[2] = 0.5*sqrt(std::max(1-R(0,0)-R(1,1)+R(2,2),0.0))*(R(1,0) - R(0,1) >= 0 ? 1 : -1);
  q[3] = 0.5*sqrt(std::max(1+R(0,0)+R(1,1)+R(2,2),0.0));

  return q;
}

inline Matrix<3,1> errFromRot(const Matrix<3,3>& R) {
  Matrix<3,1> q;
  q(0,0) = R(2,1) - R(1,2);
  q(1,0) = R(0,2) - R(2,0);
  q(2,0) = R(1,0) - R(0,1);

  double r = sqrt(q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0));
  double t = R(0,0) + R(1,1) + R(2,2) - 1;

  if (r == 0) {
    return zeros<3,1>();
  } else {
    return q * (atan2(r, t) / r);
  }
}

inline Matrix<3,3> rotFromQuat(const Matrix<4,1>& q) {
  double w = q[3];
  double x = q[0];
  double y = q[1];
  double z = q[2];

  double Nq = w*w + x*x + y*y + z*z;
  double s;
  if (Nq > 0.0) {
    s = 2/Nq;
  } else {
    s = 0.0;
  }
  double X = x*s; double Y = y*s; double Z = z*s;
  double wX = w*X; double wY = w*Y; double wZ = w*Z;
  double xX = x*X; double xY = x*Y; double xZ = x*Z;
  double yY = y*Y; double yZ = y*Z; double zZ = z*Z;

  Matrix<3,3> R;
  R(0,0) = 1.0-(yY+zZ); R(0,1) = xY-wZ;       R(0,2) = xZ+wY;
  R(1,0) = xY+wZ;       R(1,1) = 1.0-(xX+zZ); R(1,2) = yZ-wX;
  R(2,0) = xZ-wY;       R(2,1) = yZ+wX;       R(2,2) = 1.0-(xX+yY);
    
  return R;
}

inline Matrix<3,3> skewSymmetric(const Matrix<3>& vector) {
  Matrix<3,3> result = zeros<3,3>();
  result(0,1) = -vector[2]; result(0,2) = vector[1];
  result(1,0) = vector[2];    result(1,2) = -vector[0];
  result(2,0) = -vector[1]; result(2,1) = vector[0];

  return result;
}

inline double hypot(const Matrix<3>& vector) {
  double X = abs(vector[0]);
  double Y = abs(vector[1]);
  double Z = abs(vector[2]);
  if (X > Y && X > Z) {
    return X*sqrt(1.0 + (Y/X)*(Y/X) + (Z/X)*(Z/X));
  } else if (Y > Z) {
    return Y*sqrt(1.0 + (X/Y)*(X/Y) + (Z/Y)*(Z/Y));
  } else {
    return Z*sqrt(1.0 + (X/Z)*(X/Z) + (Y/Z)*(Y/Z));
  }
}


// TODO: reference additional headers your program requires here
