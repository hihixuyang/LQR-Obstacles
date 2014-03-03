// QuadrotorHoverController.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define X_DIM 16         // State dimension
#define U_DIM 4          // Control input dimension
#define Z_DIM 6          // Observation dimension

typedef Matrix<X_DIM> State;
typedef Matrix<U_DIM> Input;
typedef Matrix<Z_DIM> Observation;
typedef Matrix<3,3>   Rotation;

// world constants
double dt;              // time step, s
double gravity;         // gravity,   m/s^2
Matrix<3,1> north;      // unit vector pointing north in inertial frame
Matrix<3,1> eX, eY, eZ; // unit vectors

// quadrotor constants
double mass;            // mass, kg
Matrix<3,3> inertia;    // moment of inertia matrix 
double momentConst;     // ratio between force and moment of rotor
double dragConst;       // ratio between speed and opposite drag force
double thrust_latency;  // latency of control input, s^-1
double length;          // distance between center and rotor, m
double minForce;        // minimum force of rotor, N
double maxForce;        // maximum force of rotor, N

// derivative constants
Matrix<3,3> invInertia; // inverse of intertia matrix
double nominalInput;    // force of rotors at hover
double jStep;           // Step size for numerical derivative

// visualization
int cal_quadrotor;
int cal_scene;
int cal_floor;
int cal_ellipse;

void setupVisualization() {
  // World parameters
  dt          = 1.0/30.0;          // time step, s
  gravity     = 9.80665;           // gravity,   m/s^2
  north[0] = 0.70710678118654752440084436210485; north[1] = 0; north[2] = -0.70710678118654752440084436210485;
  eX[0] = 1; eX[1] = 0; eX[2] = 0;
  eY[0] = 0; eY[1] = 1; eY[2] = 0;
  eZ[0] = 0; eZ[1] = 0; eZ[2] = 1;
  
  // Quadrotor parameters
  mass        = 0.500;             // mass, kg  (source: paper)
  inertia     = 0.1*identity<3>(); // moment of inertia matrix 
  momentConst = 1.5e-9 / 6.11e-8;  // ratio between force and moment of rotor
  dragConst   = 0.25;
  thrust_latency = 40.0;           // latency of control input, s^-1
  length      = 0.3429/2;          // distance between center and rotor, m
  minForce    = 0.09;              // minimum force of rotor, N
  maxForce    = 3.71;              // maximum force of rotor, N
  
  // Derivative constants
  invInertia   = !inertia;			// !matrix overloaded to give the inverse of a matrix
  nominalInput = gravity*mass/4;
  jStep        = 0.0009765625;
  
  // Visualization parameters
  double beamWidth     = 0.015; // m
  double beamHeight    = 0.0065; // m
  double beamRadius    = 0.02; // m
  double motorRadius   = 0.015; // m
  double motorHeight   = 0.02; // m
  double rotorRadius   = 0.10; // m
  double rotorHeight   = 0.005; // m
  double centerSide    = 0.0889; // m
  double centerHeight  = 0.0365; // m
  double centerTopSide = 0.03; // m
  double flagLength    = 0.0508; // m
  double tileSize      = 1;  // m

  double volumeWidth   = 6.8834; // m
  double volumeLength  = 4.572;  // m
  double volumeHeight  = 3.3528;  // m

  // visualization
  CAL_Initialisation(true, true, true);
  int obj;

  // Quadrotor
  CAL_CreateGroup(&cal_quadrotor, 0, false, "QuadRotor");
  CAL_SetGroupColor(cal_quadrotor, 0.05, 0.05, 0.05);
  CAL_CreateBox(cal_quadrotor, 2*length, beamWidth, beamHeight, 0, 0, 0);
  CAL_CreateBox(cal_quadrotor, beamWidth, 2*length, beamHeight, 0, 0, 0);
  CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, length, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, -length, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, length, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, -length, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, length, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, -length, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, length, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, -length, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, -length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, -length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateBox(cal_quadrotor, centerSide, centerSide, beamHeight, 0, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, 0, 0, (float) (M_PI*0.25));
  CAL_CreateBox(cal_quadrotor, flagLength, beamWidth + 0.001, beamHeight + 0.001, length / 1.65, 0, 0, &obj);
  CAL_SetObjectColor(obj, 1, 0.15, 0);

  float flagTriangle[18] = {length / 1.65 - flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length / 1.65 - flagLength / 2, 0, -beamHeight / 2};
  CAL_CreateTriangles(cal_quadrotor, 2, flagTriangle, &obj);
  CAL_SetObjectColor(obj, 1, 0.15, 0);

  float polygon1[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
                       -sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight,
                        sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
                        sqrt(2.0)*centerSide/2, 0, 0,
                        sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
                       -sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight};
  CAL_CreatePolygon(cal_quadrotor, 6, polygon1, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon2[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
                        -sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight,
                         sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
                         sqrt(2.0)*centerSide/2, 0, 0,
                         sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
                        -sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight};
  CAL_CreatePolygon(cal_quadrotor, 6, polygon2, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon3[18] = {0, -sqrt(2.0)*centerSide/2, 0,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight,
                        0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
                        0, sqrt(2.0)*centerSide/2, 0,
                        0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight};
  CAL_CreatePolygon(cal_quadrotor, 6, polygon3, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon4[18] = {0, -sqrt(2.0)*centerSide/2, 0,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight,
                        0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
                        0, sqrt(2.0)*centerSide/2, 0,
                        0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight};
  CAL_CreatePolygon(cal_quadrotor, 6, polygon4, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  
  // Floor
  CAL_CreateGroup(&cal_floor, 0, false, "Floor");
  CAL_LoadTexture(0, "floor.png", 1);
  CAL_CreateBox(cal_floor, volumeWidth, volumeLength, 0.1, 0, 0, -0.05, &obj);
  CAL_SetObjectTexture (obj, 0, volumeWidth/tileSize, volumeLength/tileSize);

  // Ellipse
  CAL_CreateGroup (&cal_ellipse, 0, false, "Ellipse");
  CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.5);
  CAL_CreateSphere(cal_ellipse, 1, 0, 0, 0);
  
  // Volume
  CAL_CreateGroup(&cal_scene, 0, false, "Scene");
  /*CAL_CreateBox(cal_scene, volumeWidth, volumeLength, volumeHeight, 0, 0, volumeHeight / 2);
  CAL_SetGroupColor(cal_scene, 1, 1, 1, 0.2);*/
}

// Generate uniform random number 0 <= x < 1
inline double random() {
  return (double) (rand()*(RAND_MAX+1) + rand()) / (RAND_MAX*(RAND_MAX + 2));
}

// Generate random number from standard normal distribution N(0,1)
inline double normal() {
  double u, v, s(0);

  while (s == 0 || s > 1) {
    u = 2*random()-1;
    v = 2*random()-1;
    s = u*u + v*v;
  }

  return u * sqrt(-2*log(s)/s);
}

// Generate random vector from normal distribution N(mean, var)
template <size_t size>
inline Matrix<size> sampleGaussian(const Matrix<size>& mean, const Matrix<size, size>& var) {
  Matrix<size> sample;
  for (int j = 0; j < size; ++j) {
    sample[j] = normal();
  }
  Matrix<size, size> SVec, SVal;
  jacobi(var, SVec, SVal);
  for (int i = 0; i < size; ++i) {
    SVal(i,i) = sqrt(SVal(i,i));
  }
  return SVec * SVal * sample + mean;
}

// Dynamics model \dot{x} = f(x,u)
inline State f(const State& x, const Rotation& R, const Input& u) {
  State xdot;

  Matrix<3> p = x.subMatrix<3,1>(0,0);
  Matrix<3> v = x.subMatrix<3,1>(3,0);
  Matrix<3> r = x.subMatrix<3,1>(6,0);
  Matrix<3> w = x.subMatrix<3,1>(9,0);
  Matrix<4> F = x.subMatrix<4,1>(12,0);

  // \dot{p} = v
  xdot.insert(0, 0, v);

  // \dot{v} = [0,0,-g]^T + R*exp([r])*[0,0,(f_1 + f_2 + f_3 + f_4) / m]^T; 
  xdot.insert(3, 0, -gravity*eZ + R*exp(skewSymmetric(r))*((F[0]+F[1]+F[2]+F[3])/mass)*eZ); 

  // \dot{r} = w + 0.5*skewSymmetric(r)*w + (1.0/tr(~r*r))*(1.0 - 0.5*sqrt(tr(~r*r))/tan(0.5*sqrt(tr(~r*r))))*skewSymmetric(r)*(skewSymmetric(r)*w)
  double l = hypot(x.subMatrix<3,1>(6,0));
  if (0.5*l > 0.0) {
    xdot.insert(6, 0, w + 0.5*skewSymmetric(r)*w + (1.0 - 0.5*l/tan(0.5*l))*skewSymmetric(r / l)*(skewSymmetric(r / l)*w));
  } else {
    xdot.insert(6, 0, w + 0.5*skewSymmetric(r)*w);
  }
  // \dot{w} = J^{-1}*([l*(f_2 - f_4), l*(f_3 - f_1), (f_1 - f_2 + f_3 - f_4)*k_M]^T - [w]*J*w)
  xdot.insert(9, 0, invInertia*( length*(F[1] - F[3])*eX + length*(F[2] - F[0])*eY + (F[0] - F[1] + F[2] - F[3])*momentConst*eZ - skewSymmetric(w)*inertia*w));

  // \dot{f} = latency*(f^* - f)
  xdot.insert(12,0, thrust_latency*(u - F));

  return xdot;
}

inline Observation h(const State& x, const Rotation& R) {
  Observation z;

  // rate-gyros
  z.insert(0,0,x.subMatrix<3,1>(9,0));

  // tracking
  z.insert(3,0,x.subMatrix<3,1>(0,0));

  // accelerometer
  /*z[0] = -(dragConst/mass)*((R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*x[3] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*x[4] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*x[5]);
  z[1] = -(dragConst/mass)*((R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*x[3] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*x[4] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*x[5]);
  z[2] = (forceConst/mass)*x[12];*/

  // magnetic compass
  /*z[3] = (R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*north[0] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*north[1] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*north[2];
  z[4] = (R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*north[0] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*north[1] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*north[2];
  z[5] = (R(0,2) + R(0,0)*x[7] - R(0,1)*x[6])*north[0] + (R(1,2) + R(1,0)*x[7] - R(1,1)*x[6])*north[1] + (R(2,2) + R(2,0)*x[7] - R(2,1)*x[6])*north[2];*/

  return z;
}

inline Matrix<X_DIM, X_DIM> Jacobian_fx(const State& x, const Rotation& R, const Input& u) {
  Matrix<X_DIM,X_DIM> F;
  Matrix<X_DIM> xr(x), xl(x);
  for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += jStep; xl[i] -= jStep;
		F.insert(0,i, (f(xr, R, u) - f(xl, R, u)) / (2*jStep));
    xr[i] = xl[i] = x[i];
	}
  return F;
}

inline Matrix<X_DIM, U_DIM> Jacobian_fu(const State& x, const Rotation& R, const Input& u) {
  Matrix<X_DIM,U_DIM> G;
  Matrix<U_DIM> ur(u), ul(u);
  for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += jStep; ul[i] -= jStep;
		G.insert(0,i, (f(x, R, ur) - f(x, R, ul)) / (2*jStep));
    ur[i] = ul[i] = u[i];
	}
  return G;
}

inline Matrix<Z_DIM, X_DIM> Jacobian_hx(const State& x, const Rotation& R) {
  Matrix<Z_DIM, X_DIM> H;
  Matrix<X_DIM> xr(x), xl(x);
  for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += jStep; xl[i] -= jStep;
		H.insert(0,i, (h(xr, R) - h(xl, R)) / (2*jStep));
    xr[i] = xl[i] = x[i];
	}
  return H;
}

void visualize(const State& xTrue, const Rotation& RTrue, const State& xEst, const Rotation& REst, const Matrix<X_DIM,X_DIM>& P, double t) {
  float p[3] = {(float) xTrue[0], (float) xTrue[1], (float) xTrue[2]};
  Matrix<4> q = quatFromRot(RTrue);
  float o[4] = {(float) q[0], (float) q[1], (float) q[2], (float) q[3]};
  CAL_AddGroupKeyState(cal_quadrotor, (float) t, p, o);

  Matrix<3,3> V, E;
  jacobi(P.subMatrix<3,3>(0,0), V, E);
  q = quatFromRot(V);
  float pos[3] = {(float) xEst[0], (float) xEst[1], (float) xEst[2]};
  float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
  float scale[3] = {(float) (2*sqrt(E(0,0))), (float) (2*sqrt(E(1,1))), (float) (2*sqrt(E(2,2)))};
  CAL_AddGroupKeyState(cal_ellipse, (float) t, pos, quat, scale);
}

inline void linearizeDiscretize(const State& x, const Rotation& R, const Input& u, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,1>& c) {
  // Linearize:  \dot{x} = f(\hat{x}, R, \hat{u}) + F(x - \hat{x}) + G(u - \hat{u})
  // Discretize: (x - \hat{x})_{t+1} = A*(x - \hat{x})_t + B*(u - \hat{u})_t + c

  Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
  Matrix<X_DIM,U_DIM> G = Jacobian_fu(x,R,u);
  Matrix<X_DIM,1>  xDot = f(x,R,u);
  
  A = exp(dt*F);
 
  // Simpson's integration approximation of Int[0,dt] exp(F*T) dT
  Matrix<X_DIM,X_DIM> Int = (dt/6.0) * (identity<X_DIM>() + 4.0*exp(0.5*dt*F) + A);
  
  B = Int * G;
  c = Int * xDot;
}

inline void propagate(State& x, Rotation& R, const Input& u, const Matrix<X_DIM,X_DIM>& M) {
  Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
  State xDot = f(x,R,u);
  
  Matrix<X_DIM,X_DIM> A = exp(dt*F);
  Matrix<X_DIM,X_DIM> A2 = exp((dt*0.5)*F);
  Matrix<X_DIM,X_DIM> MM = (dt/6) * (M + 4*A2*M*~A2 + A*M*~A);

  x = x + (dt/6)*(xDot + 4*(A2*xDot) + A*xDot) + sampleGaussian(zeros<X_DIM,1>(), MM);
  
  // reset rotation-error in xNew into RNew
  R = R*exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
  x[6] = 0; x[7] = 0; x[8] = 0;
}

inline void kalmanFilter1(State& x, Rotation& R, const Input& u, const Matrix<X_DIM,X_DIM>& M, Matrix<X_DIM,X_DIM>& P) {
  Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
  State xDot = f(x,R,u);
  
  Matrix<X_DIM,X_DIM> A = exp(dt*F);
  Matrix<X_DIM,X_DIM> A2 = exp((dt*0.5)*F);

  // update estimate of state
  x = x + (dt/6)*(xDot + 4*(A2*xDot) + A*xDot);
  
  // update state covariance
  Matrix<X_DIM,X_DIM> MM = (dt / 6) * (M + 4*A2*M*~A2 + A*M*~A);
  P = A*P*~A + MM;

  // reset rotation-error in x into R
  R = R*exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
  x[6] = 0; x[7] = 0; x[8] = 0;
}

inline void kalmanFilter2(State& x, Rotation& R, const Observation& z, const Matrix<Z_DIM, Z_DIM>& N, Matrix<X_DIM, X_DIM>& P) {
  // incorporate measurement
  Matrix<Z_DIM,X_DIM> H = Jacobian_hx(x,R);
  Matrix<X_DIM,Z_DIM> K = P*~H*!(H*P*~H + N);
  
  x = x + K*(z - h(x,R));
  P = (identity<X_DIM>() - K*H)*P;

  // reset rotation-error in x into R
  R = R*exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
  x[6] = 0; x[7] = 0; x[8] = 0;
}

inline void controlMatrices(const Input& uGoal, const Matrix<X_DIM, X_DIM>& Q, const Matrix<U_DIM, U_DIM>& R, Matrix<U_DIM,X_DIM>& L, Matrix<U_DIM,X_DIM>& E, Matrix<U_DIM,1>& l) {
  Matrix<X_DIM,X_DIM> A;
  Matrix<X_DIM,U_DIM> B;
  Matrix<X_DIM,1> c;
  
  // Linearize about steady state: velocity, roll, pitch, angular velocity, control are 0. Position and yaw are free.
  State xHat;
  Rotation RHat;

  xHat[0] = 0;     xHat[1] = 0;     xHat[2] = 0;
  xHat[3] = 0;     xHat[4] = 0;     xHat[5] = 0;
  xHat[6] = 0;     xHat[7] = 0;     xHat[8] = 0;
  xHat[9] = 0;     xHat[10] = 0;    xHat[11] = 0;
  xHat[12] = uGoal[0]; xHat[13] = uGoal[1]; xHat[14] = uGoal[2]; xHat[15] = uGoal[3];

  RHat = identity<3>();

  linearizeDiscretize(xHat, RHat, uGoal, A, B, c);

  // Solve riccati equation to find S (naive implementation)
  Matrix<X_DIM, X_DIM> S = Q;
  for (int i = 0; i < 200; ++i) {
    S = Q + ~A*S*A - ~A*S*B*!(R + ~B*S*B)*~B*S*A;
  }

  Matrix<X_DIM,X_DIM> T = !(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * Q;
  Matrix<X_DIM,1>     v = !(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * (~A*S*B*!(R + ~B*S*B)*~B*S - ~A*S)*c;
  
  // feedback law: u = L*xTilde + E*xTildeGoal + l
  L = -!(R + ~B*S*B)*~B*S*A;
  E = -!(R + ~B*S*B)*~B*T;
  l = -!(R + ~B*S*B)*(~B*S*c + ~B*v); // should be zero for steady-state linearization
}

inline Input riccatiControllerSteady(const State& x, const Rotation& R0, const State& xGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<U_DIM,X_DIM>& L, const Matrix<U_DIM,X_DIM>& E, const Matrix<U_DIM,1>& l) {
  
  // find rotation such that yaw = 0 (express goal and current state relative to linearization state)
  Rotation RLocal;
  
  Matrix<3> zVecR0, zVec;
  zVecR0[0] = R0(0,2); zVecR0[1] = R0(1,2); zVecR0[2] = R0(2,2);
  zVec[0]   = 0;       zVec[1]   = 0;       zVec[2]   = 1;

  Matrix<3> axis = skewSymmetric(zVecR0)*zVec;
  double sinangle = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
  double angle = asin(sinangle);
  if (sinangle != 0) {
    axis = axis * (angle / sinangle);
  }

  RLocal = exp(skewSymmetric(axis)) * R0;
  
  State xTildeGoal; //= xGoal - xHat;
  xTildeGoal.insert(0,0, ~RLocal * xGoal.subMatrix<3,1>(0,0));
  xTildeGoal.insert(3,0, ~RLocal * xGoal.subMatrix<3,1>(3,0));
  xTildeGoal.insert(6,0, errFromRot(~RLocal*RGoal));
  xTildeGoal.insert(9,0, xGoal.subMatrix<3,1>(9,0));
  xTildeGoal.insert(12,0, xGoal.subMatrix<4,1>(12,0) - uGoal);

  State xTilde; //= x - xHat;
  xTilde.insert(0,0, ~RLocal * x.subMatrix<3,1>(0,0));
  xTilde.insert(3,0, ~RLocal * x.subMatrix<3,1>(3,0));
  xTilde.insert(6,0, errFromRot(~RLocal*R0)); //xTilde.insert(6,0, -(~R0*axis));
  xTilde.insert(9,0, x.subMatrix<3,1>(9,0));
  xTilde.insert(12,0, x.subMatrix<4,1>(12,0) - uGoal);
  
  return uGoal + L*xTilde + E*xTildeGoal + l;
}

inline Input riccatiControllerCurrent(const State& x, const Rotation& R0, const State& xGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<X_DIM, X_DIM>& Q, const Matrix<U_DIM, U_DIM>& R, Matrix<X_DIM, X_DIM>& S) {
  Matrix<X_DIM,X_DIM> A;
  Matrix<X_DIM,U_DIM> B;
  Matrix<X_DIM,1> c;
  
  // Linearize about current state.
  State xHat = x;
  Rotation RHat = R0;

  linearizeDiscretize(xHat, RHat, uGoal, A, B, c);

  // Solve riccati equation to find S (naive implementation)
  for (int i = 0; i < 100; ++i) {
    S = Q + ~A*(S - S*B*!(R + ~B*S*B)*~B*S)*A;
  }
  
  Matrix<X_DIM,X_DIM> T = !(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * Q;
  Matrix<X_DIM,1>     v = -!(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * ~A*(S - S*B*!(R + ~B*S*B)*~B*S)*c;
  
  // feedback law: u = L*xTilde + E*xTildeGoal + l
  Matrix<U_DIM,X_DIM> L = -!(R + ~B*S*B)*~B*S*A;
  Matrix<U_DIM,X_DIM> E = -!(R + ~B*S*B)*~B*T;
  Matrix<U_DIM,1>     l = -!(R + ~B*S*B)*~B*(S*c + v);

  // solve feedforward vector given S
  State xTildeGoal = xGoal - xHat;
  xTildeGoal.insert(6,0, errFromRot(~RHat*RGoal));

  State xTilde = zeros<X_DIM,1>();
  
  return uGoal + L*xTilde + E*xTildeGoal + l;
}

int _tmain(int argc, _TCHAR* argv[])
{
  setupVisualization();

  State xInit, xGoal;
  Input uGoal;

  xInit[0] = 0; xInit[1] = 0; xInit[2] = 1;						// pos
  xInit[3] = 0;  xInit[4] = 0;  xInit[5] = 0;					// vel
  xInit[6] = 0;  xInit[7] = 0;  xInit[8] = 0;					// or relative to Rot
  xInit[9] = 0; xInit[10] = 0; xInit[11] = 0;					// ang speed
  xInit[12] = xInit[13] = xInit[14] = xInit[15] = nominalInput; // thrust

  Rotation RotInit = identity<3>();
  
  xGoal[0] = 2; xGoal[1] = 3; xGoal[2] = 5;						// pos
  xGoal[3] = 0;  xGoal[4] = 0;  xGoal[5] = 5;					// vel
  xGoal[6] = 0;  xGoal[7] = 0;  xGoal[8] = 0;					// rotation vector
  xGoal[9] = 0; xGoal[10] = 0; xGoal[11] = 0;					// ang speed
  xGoal[12] = xGoal[13] = xGoal[14] = xGoal[15] = nominalInput;	// thrust

  Rotation RotGoal = identity<3>();
  RotGoal(0,0) = 1; RotGoal(0,1) = 0;
  RotGoal(1,0) = 0; RotGoal(1,1) = 1;

  uGoal[0] = uGoal[1] = uGoal[2] = uGoal[3] = nominalInput;

  Matrix<X_DIM, X_DIM> Q = zeros<X_DIM,X_DIM>();
  Q(0,0) = Q(1,1) = Q(2,2) = 10;
  //Q(0,0) = Q(1,1) = 0; Q(2,2) = 0;
  Q(3,3) = Q(4,4) = Q(5,5) = 1;
  //Q(3,3) = Q(4,4) = Q(5,5) = 10;
  Q(6,6) = Q(7,7) = 1; Q(8,8) = 20;
  Q(9,9) = 1; Q(10,10) = 1; Q(11,11) = 1;
  
  Matrix<U_DIM, U_DIM> R = 20*identity<U_DIM>();
  
  Matrix<X_DIM,X_DIM> P = 0.0001 * 0.0001 * identity<X_DIM>(); // initial state variance
  P(0,0) = P(1,1) = 0.1*0.1; P(2,2) = 0.0001*0.0001; // position (z very certain, because on floor)
  P(3,3) = P(4,4) = P(5,5) = 0.0001*0.0001; // velocity (very certain, because on floor)
  P(6,6) = P(7,7) = P(8,8) = 0.0001*0.0001; // orientation (zero by definition)
  P(9,9) = P(10,10) = P(11,11) = 0.0001*0.0001; // angular speed (very certain, because on floor)
  P(12,12) = P(13,13) = P(14,14) = P(15,15) = 0.1*0.1; // force

  Matrix<X_DIM,X_DIM> M = 0.0001 * identity<X_DIM>(); // motion noise variance
  M(0,0) = M(1,1) = M(2,2) = 0.0001*0.0001; // position
  M(3,3) = M(4,4) = M(5,5) = 0.01*0.01; // velocity
  M(6,6) = M(7,7) = M(8,8) = 0.0001*0.0001; // orientation
  M(9,9) = M(10,10) = M(11,11) = 0.01*0.01; // angular speed
  M(12,12) = P(13,13) = P(14,14) = P(15,15) = 0.01*0.01; // force

  Matrix<Z_DIM,Z_DIM> N = identity<Z_DIM>(); // observation noise variance
  N(0,0) = N(1,1) = N(2,2) = 0.001*0.001; // gyro
  N(3,3) = N(4,4) = N(5,5) = 0.01*0.01; // position (in units)
    
  clock_t beginTime = clock();

  State x = xInit;
  Rotation Rot = RotInit;

  State xTrue = sampleGaussian(xInit, P);
  Rotation RotTrue = RotInit;
  RotTrue = RotTrue*exp(skewSymmetric(xTrue.subMatrix<3,1>(6,0)));
  xTrue[6] = 0; xTrue[7] = 0; xTrue[8] = 0;

  Matrix<U_DIM,X_DIM> L; Matrix<U_DIM,X_DIM> E; Matrix<U_DIM,1> l;
  controlMatrices(uGoal, Q, R, L, E, l);

  Matrix<X_DIM,X_DIM> S = Q;

  for (size_t t = 0; t < 320; ++t) {
    Input u = riccatiControllerSteady(x, Rot, xGoal, RotGoal, uGoal, L, E, l);
    //Input u = riccatiControllerCurrent(x, Rot, xGoal, RotGoal, uGoal, Q, R, S);

    // Clamp control input to minimum and maximum value 
    /*for (size_t i = 0; i < U_DIM; ++i) {
      if (u[i] < minForce) {
        std::cout << "_" << u[i] << " ";
        u[i] = minForce;    
      }
      if (u[i] > maxForce) {
        std::cout << "^" << u[i] << " ";
        u[i] = maxForce;
      }
    }*/

    // Send control input to quadrotor
    propagate(xTrue, RotTrue, u, M);
    // Process control input in time update of Kalman filter
    kalmanFilter1(x, Rot, u, M, P);

    // Receive observation
    Observation z = sampleGaussian(h(xTrue, RotTrue), N);
    // Process observation in measurement of Kalman filter
    kalmanFilter2(x, Rot, z, N, P);

    visualize(xTrue, RotTrue, x, Rot, P, t*dt);
  }

  std::cout << (clock() - beginTime) / (double) (320*CLOCKS_PER_SEC) << std::endl;

  int k;
  std::cin >> k;
  
  CAL_End();
 	return 0;
}