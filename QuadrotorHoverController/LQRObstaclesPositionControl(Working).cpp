/*
Daman Bareiss
Algorithmic Robotics Laboratory
University of Utah
(405)-642-9754
daman.bareiss@gmail.com
*/


#include "simulator2.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "gjk.h"
#include "Vector3.h"
#include "Windows.h"

/*
NUM_QUADS		- The number of quadrotors that will be in the simulation. Select a number > 1. 
				      Currently must be 6 or less but on line 1172 you can see how to add the initial condition for more than 6
FREQUENCY		- Set the frequency of the simulation. 30Hz has been sufficient in all my experiments
TIME_STEPS		- The length the simulation will run in seconds
OBSTACLE_STEPS	- The lengtho of the time horizon in seconds
NUM_POINTS		- The number of points that are sampled around the quadrotor's ellipse and copied to build the LQR-Obstacle
XYRADIUS		- The minor radii of the quadrotors' bounding ellipsoid
ZRADIUS			- The major radii of the quadrotors' bounding ellipsoid
*/

#define NUM_QUADS 6		
#define FREQUENCY 30
#define TIME_STEPS 7*FREQUENCY
#define OBSTACLE_STEPS 30	
#define NUM_POINTS 50		
#define XYRADIUS 0.26
#define ZRADIUS 0.75 
bool drawEllipses = false;

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
double length_;          // distance between center and rotor, m
double minForce;        // minimum force of rotor, N
double maxForce;        // maximum force of rotor, N

// derivative constants
Matrix<3,3> invInertia; // inverse of intertia matrix
double nominalInput;    // force of rotors at hover
double jStep;           // Step size for numerical derivative

// visualization
int rotordummy;
int ellipsedummy;
int cal_scene;
int cal_floor;
int cal_ellipse;

int counter=0;

// LQR and Noise Matrices
Matrix<X_DIM,X_DIM> Qx;
Matrix<V_DIM,V_DIM> Qv;
Matrix<V_DIM,V_DIM> Qp;
Matrix<U_DIM,U_DIM> R;
Matrix<U_DIM,U_DIM> Rp;
Matrix<X_DIM,X_DIM> M;
Matrix<Z_DIM,Z_DIM> N;

// Quadrotor Class
class Quadrotor {
public:
  State xGoal;		// Desired State
  Position pGoal;	// Desired Position
  Velocity vGoal;	// Desired Velocity (from Position Controller)
  Velocity newV;	// New, collision free velocity (from Linear Programming)
  Input uGoal;		// Goal low-level control input from Velocity LQR control
  Rotation RotGoal;
  
  int cal_quadrotor; // For visualization
  int cal_ellipse;	 // For visualization

  State x;			// State with noise
  Rotation Rot;		// Rotation with noise
  State xTrue;		// State without noise
  Rotation RotTrue;	// Rotation without noise
  Matrix<U_DIM,X_DIM> L;	// Used in Velocity-LQR
  Matrix<U_DIM,V_DIM> E;	// Used in Velocity-LQR
  Matrix<U_DIM,1> l;		// ell, used in Velocity-LQR
  Matrix<V_DIM,X_DIM> Lh;	// Used in Position-LQR
  Matrix<V_DIM,V_DIM> Eh;	// Used in Position-LQR

  Matrix<X_DIM,X_DIM> P;	// Noise covariance on state

  std::vector< Matrix<3,1> > ellipsoids;		// Vector of points making up the whole obstacle
  std::vector< Matrix<3,1> > reachablePoints;	// Vector of points from the obstacle that are reachable (within control input constraints)
  Matrix<3,3> Transform;						// Transformation matrix PG^-1
  Matrix<3,1> Translate;						// Translation Matrix
  std::vector< Matrix<4,1> > planes;			// Convex hull facets' plane equation coefficients
  std::vector< Matrix<3,1> > normals;			// Convex hull normal vector
  std::vector< Matrix<3,3> > vertex;			// Convex hull facet's plane vertex
  Matrix<4,1> shortestDistance;					// Escape vector. Elements 0-2 = direction, 3 = magnitude 

  // Initialize a quadrotor with the given initial conditions
  void setupQuadrotors(const State& xinit, const Rotation& rotinit, const Matrix<X_DIM,X_DIM> Pinit, const State& xfinal, const Position& pfinal, const Input& uHover, const Rotation& rotfinal) {
    x = xinit;
    Rot = rotinit;
    xGoal = xfinal;
    pGoal = pfinal;
	vGoal = zeros<V_DIM,1>();
    uGoal = uHover;
    RotGoal = rotfinal;
    P = Pinit;

	ellipsoids.reserve(NUM_POINTS*OBSTACLE_STEPS);
	reachablePoints.reserve(NUM_POINTS*OBSTACLE_STEPS);

    xTrue = sampleGaussian(x,P);
    RotTrue = Rot;
    RotTrue = RotTrue * exp(skewSymmetric((xTrue).subMatrix<3,1>(6,0)));
    xTrue[6] = 0; xTrue[7] = 0; xTrue[8] = 0;
  }
  // Creates the visualization of the quadrotor
  void setupQuadVisualization() {
    CAL_CloneGroup(&cal_quadrotor, rotordummy, 0, false, "Quadrotor");
    CAL_CloneGroup(&cal_ellipse, ellipsedummy, 0, false, "ellipse");
  }
  // Sets the visualization for a given time step t
  void visualize(double t) {
    float p[3] = {(float) xTrue[0], (float) xTrue[1], (float) xTrue[2]};
    Matrix<4> q = quatFromRot(RotTrue);
    float o[4] = {(float) q[0], (float) q[1], (float) q[2], (float) q[3]};
    CAL_AddGroupKeyState(cal_quadrotor, (float) t, p, o);
	CAL_AddGroupKeyState(cal_ellipse, (float) t, p);
 
    /*Matrix<3,3> V, E;
    jacobi(P.subMatrix<3,3>(0,0), V, E);
    q = quatFromRot(V);
    float pos[3] = {(float) x[0], (float) x[1], (float) x[2]};
    float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
    float scale[3] = {(float) (2*sqrt(E(0,0))), (float) (2*sqrt(E(1,1))), (float) (2*sqrt(E(2,2)))};
    CAL_AddGroupKeyState(cal_ellipse, (float) t, pos, quat, scale);
	CAL_SetGroupPosition(cal_ellipse, p[0], p[1], p[2]);
    CAL_SetGroupQuaternion(cal_ellipse, q[0], q[1], q[2], q[3]);*/
  }

  // Finds the linearized-discretized state and control matrices
  void findMatrices(Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,1>& c) {
    controlMatrices(uGoal, xGoal, A, B, c, L, E, l, Lh, Eh);
  }
  
  // Find the low level input through the Velocity-LQR
  Input findU() {
    return riccatiControllerSteady(x, Rot, vGoal, RotGoal, uGoal, L, E, l);
  }

  // Find the goal velocity through the Position-LQR
   Velocity findVGoal() {
    return riccatiControllerSteadyPosition(x,Rot, pGoal, RotGoal, uGoal, Lh, Eh);
  }

  // Send the input to the robot
  void propagateU( const Input& u ) {
    propagate(xTrue, RotTrue, u, M);
  }
};

void setup() {
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
  length_      = 0.3429/2;          // distance between center and rotor, m
  minForce    = 0.09;              // minimum force of rotor, N
  maxForce    = 3.71;              // maximum force of rotor, N
  
  // Derivative constants
  invInertia   = !inertia;
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

  // visualization setup
  CAL_Initialisation(true, true, true);
  CAL_SetViewParams(0,0,0,20,0,0,0, 0, -1, 0);
  int obj;
  
  // Quadrotor
  CAL_CreateGroup(&rotordummy, 0, false, "QuadRotor");
  CAL_SetGroupColor(rotordummy, 0.05, 0.05, 0.05);
  CAL_CreateBox(rotordummy, 2*length_, beamWidth, beamHeight, 0, 0, 0);
  CAL_CreateBox(rotordummy, beamWidth, 2*length_, beamHeight, 0, 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, length_, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, -length_, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, 0, length_, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, 0, -length_, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, length_, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, -length_, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, 0, length_, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, 0, -length_, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, length_, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, -length_, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, 0, length_, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, 0, -length_, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateBox(rotordummy, centerSide, centerSide, beamHeight, 0, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, 0, 0, (float) (M_PI*0.25));
  CAL_CreateBox(rotordummy, flagLength, beamWidth + 0.001, beamHeight + 0.001, length_ / 1.65, 0, 0, &obj);
  CAL_SetObjectColor(obj, 1, 0.15, 0);

  float flagTriangle[18] = {length_ / 1.65 - flagLength / 2, 0, -beamHeight / 2,
                            length_ / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length_ / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length_ / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length_ / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length_ / 1.65 - flagLength / 2, 0, -beamHeight / 2};
  CAL_CreateTriangles(rotordummy, 2, flagTriangle, &obj);
  //CAL_CreateTriangles(rotordummy,2,
  CAL_SetObjectColor(obj, 1, 0.15, 0);

  float polygon1[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
                       -sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight,
                        sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
                        sqrt(2.0)*centerSide/2, 0, 0,
                        sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
                       -sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight};
  CAL_CreatePolygon(rotordummy, 6, polygon1, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon2[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
                        -sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight,
                         sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
                         sqrt(2.0)*centerSide/2, 0, 0,
                         sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
                        -sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight};
  CAL_CreatePolygon(rotordummy, 6, polygon2, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon3[18] = {0, -sqrt(2.0)*centerSide/2, 0,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight,
                        0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
                        0, sqrt(2.0)*centerSide/2, 0,
                        0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight};
  CAL_CreatePolygon(rotordummy, 6, polygon3, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
  float polygon4[18] = {0, -sqrt(2.0)*centerSide/2, 0,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight,
                        0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
                        0, sqrt(2.0)*centerSide/2, 0,
                        0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
                        0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight};
  CAL_CreatePolygon(rotordummy, 6, polygon4, &obj);
  CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);

  // Floor
  /*CAL_CreateGroup(&cal_floor, 0, false, "Floor");
  CAL_LoadTexture(0, "floor.png", 1);
  CAL_CreateBox(cal_floor, volumeWidth, volumeLength, 0.1, 0, 0, -0.05, &obj);
  CAL_SetObjectTexture (obj, 0, volumeWidth/tileSize, volumeLength/tileSize);*/

  // Ellipse
  CAL_CreateGroup (&ellipsedummy, 0, false, "Ellipse");
  CAL_SetGroupColor(ellipsedummy, 0, 1, 0, 0.5);
  CAL_CreateSphere(ellipsedummy, 1, 0, 0, 0, &obj);
  CAL_SetObjectScaling(obj, XYRADIUS, XYRADIUS, ZRADIUS);
  
  // Volume
  CAL_CreateGroup(&cal_scene, 0, false, "Scene");
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
  xdot.insert(9, 0, invInertia*( length_*(F[1] - F[3])*eX + length_*(F[2] - F[0])*eY + (F[0] - F[1] + F[2] - F[3])*momentConst*eZ - skewSymmetric(w)*inertia*w));

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
  z[2] = (forceConst/mass)*x[12];     */ 

  // magnetic compass
  /*z[3] = (R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*north[0] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*north[1] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*north[2];
  z[4] = (R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*north[0] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*north[1] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*north[2];
  z[5] = (R(0,2) + R(0,0)*x[7] - R(0,1)*x[6])*north[0] + (R(1,2) + R(1,0)*x[7] - R(1,1)*x[6])*north[1] + (R(2,2) + R(2,0)*x[7] - R(2,1)*x[6])*north[2]; */


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

// Solve for the LQR control matrices L,E,l, Lh, and Eh
inline void controlMatrices(const Input& uGoal, const Matrix<X_DIM,1>& xstar, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,1>& c, Matrix<U_DIM,X_DIM>& L, Matrix<U_DIM,V_DIM>& E, Matrix<U_DIM,1>& l, Matrix<V_DIM,X_DIM>& Lh, Matrix<V_DIM,V_DIM>& Eh) {

  Matrix<V_DIM,X_DIM> V = zeros<V_DIM,X_DIM>();
  V(0,3) = V(1,4) = V(2,5) = 1;
  Matrix<V_DIM,X_DIM> Pmap = zeros<V_DIM,X_DIM>();
  Pmap(0,0) = Pmap(1,1) = Pmap(2,2) = 1;

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
  Matrix<X_DIM, X_DIM> S = ~V*Qv*V;
  Matrix<X_DIM, V_DIM> T = -~V*Qv;
  for (int i = 0; i < 300; ++i) {
     T = -~V*Qv + ~A*T - ~A*S*B*!(R + ~B*S*B)*~B*T;
     S = ~V*Qv*V + Qx + ~A*S*A - ~A*S*B*!(R + ~B*S*B)*(~B*S*A);
  }

   Matrix<X_DIM,1>     a = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>())*(Qx*xstar - ~A*S*c + ~A*S*B*!(R+~B*S*B)*~B*S*c);

   // feedback law: uGoal = L*x + E*vGoal + l
   L = -!(R + ~B*S*B)*~B*S*A;
   E = -!(R + ~B*S*B)*~B*T;
   l = -!(R + ~B*S*B)*(~B*S*c + ~B*a);

   Matrix<X_DIM,X_DIM> Qptilde = ~Pmap*Qp*Pmap + 0.25*~L*R*L;
   Matrix<V_DIM,V_DIM> Rtilde = 0.25*~E*R*E;
   Matrix<V_DIM,X_DIM> Ptilde = 0.25*~E*R*L;
   Matrix<X_DIM,X_DIM> Atilde = A+B*L;
   Matrix<X_DIM,V_DIM> Btilde = B*E;

  // Solve riccati equation to find S (naive implementation)
  Matrix<X_DIM, X_DIM> Stilde = Qptilde;
  Matrix<X_DIM, V_DIM> Ttilde = -~Pmap*Qp;
  for (int i = 0; i < 300; ++i) {
     Ttilde = -~Pmap*Qp + ~Atilde*Ttilde - (~Ptilde + ~Atilde*Stilde*Btilde)*!(Rtilde + ~Btilde*Stilde*Btilde)*~Btilde*Ttilde; 
     Stilde = Qptilde + ~Atilde*Stilde*Atilde - (~Ptilde + ~Atilde*Stilde*Btilde)*!(Rtilde+~Btilde*Stilde*Btilde)*(Ptilde+~Btilde*Stilde*Atilde);
  }
  Matrix<X_DIM,1> atilde = zeros<X_DIM,1>();

   // feedback law: vGoal = L*x + E*pGoal
   Lh = -(!(Rtilde+~Btilde*Stilde*Btilde))*(Ptilde+~Btilde*Stilde*Atilde);
   Eh = -(!(Rtilde + ~Btilde*Stilde*Btilde))*(~Btilde*Ttilde);
}

// Find u from vGoal
inline Input riccatiControllerSteady(const State& x, const Rotation& R0, const Velocity& vGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<U_DIM,X_DIM>& L, const Matrix<U_DIM,V_DIM>& E, const Matrix<U_DIM,1>& l) {
  
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
  
  Velocity vTildeGoal; //= xGoal - xHat;
  vTildeGoal.insert(0,0, ~RLocal * vGoal);

  State xTilde; //= x - xHat;
  xTilde.insert(0,0, ~RLocal * x.subMatrix<3,1>(0,0));
  xTilde.insert(3,0, ~RLocal * x.subMatrix<3,1>(3,0));
  xTilde.insert(6,0, errFromRot(~RLocal*R0)); //xTilde.insert(6,0, -(~R0*axis));
  xTilde.insert(9,0, x.subMatrix<3,1>(9,0));
  xTilde.insert(12,0, x.subMatrix<4,1>(12,0) - uGoal);
  
  return uGoal + L*xTilde + E*vTildeGoal + l;
}

// Find goal velocity from goal position
inline Velocity riccatiControllerSteadyPosition(const State& x, const Rotation& R0, const Position& pGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<V_DIM,X_DIM>& Lh, const Matrix<V_DIM,V_DIM>& Eh) {
  
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
  
  Position pTildeGoal; //= xGoal - xHat;
  pTildeGoal.insert(0,0, ~RLocal * pGoal);

  State xTilde; //= x - xHat;
  xTilde.insert(0,0, ~RLocal * x.subMatrix<3,1>(0,0));
  xTilde.insert(3,0, ~RLocal * x.subMatrix<3,1>(3,0));
  xTilde.insert(6,0, errFromRot(~RLocal*R0)); //xTilde.insert(6,0, -(~R0*axis));
  xTilde.insert(9,0, x.subMatrix<3,1>(9,0));
  xTilde.insert(12,0, x.subMatrix<4,1>(12,0) - uGoal);
  
  return RLocal*(Lh*xTilde + Eh*pTildeGoal);
}

// GJK Functions
int
      gjk_num_g_test,     /* how many times the G-test is performed -- the
                             same as the number of main-loop iterations */
      gjk_num_simplices,  /* how many times the simplex routine
                             was called */
      gjk_num_backups,    /* how many times (if ever!) the GJK backup
                             procedure was called */
      gjk_num_dot_products, /* how many dot-product operations are called */
      gjk_num_support_dp, /* how many dot-product operations are called
			      whilst executing the support function */
      gjk_num_other_ops; /* how many other mults and divides are called */

void apply_trans(  Transform t, REAL * src, REAL * tgt)
{
  int i;

  if ( t==0 )
    for ( i=0 ; i<DIM ; i++ )
      tgt[i] = src[i];
  else {
    for ( i=0 ; i<DIM ; i++ )
      tgt[i] = t[i][DIM] + OTHER_DOT_PRODUCT( t[i], src);
  }
  return;
}
void
apply_rot_transpose( Transform t, REAL * src, REAL * tgt)
{
  int i;

  if ( t==0 )
    for ( i=0 ; i<DIM ; i++ )
      tgt[i] = src[i];
  else {
    for ( i=0 ; i<DIM ; i++ )
      tgt[i] = DO_MULTIPLY( t[0][i], src[0]) + DO_MULTIPLY( t[1][i], src[1])
	             + DO_MULTIPLY( t[2][i], src[2]);
  }
  return;
}

// Find the F and G matrices frmo dynamics
inline void findFG(const Matrix<X_DIM,X_DIM>& A, const Matrix<X_DIM,U_DIM>& B, const Matrix<U_DIM,X_DIM>& L, const Matrix<U_DIM,V_DIM>& E, Matrix<X_DIM,X_DIM>& F, Matrix<X_DIM,V_DIM>& G) {
	Matrix<X_DIM,X_DIM> Atilde;
	Atilde = A+(B*L);

	Matrix<X_DIM,V_DIM> Btilde;
	Btilde = (B*E);

	F = Atilde*F;
	G = Atilde*G + Btilde;
}

// Creates an ellipsoid of NUM_POINTS points
inline void createSpheres(std::vector< Matrix<3,1> >& points) {
	double dlong = M_PI*(3.0-sqrt(5.0));
	double dz = 2.0/NUM_POINTS;
	double longt = 0;
	double z = 1.0 -dz/2.0;
	Matrix<3,1> tempPoint;
	
	for(int i = 0; i<NUM_POINTS; i++) {
		tempPoint[0] = 2*XYRADIUS*cos(longt)*sqrt(1-z*z);
		tempPoint[1] = 2*XYRADIUS*sin(longt)*sqrt(1-z*z);
		tempPoint[2] = 2*ZRADIUS*z;
		points.push_back(tempPoint);
		z -= dz;
		longt += dlong;
	}
}

// Takes the ellipsoid of and transform it into the LQR-Obstacle for a single time step
inline void createObstacle(std::vector< Matrix<3,1> >& ellipsoids, const std::vector< Matrix<3,1> >& points, Matrix<3,3>& Transform, Matrix<3,1>& Translate, const size_t& t, const Matrix<X_DIM,X_DIM>& F, const Matrix<X_DIM,V_DIM>& G, const Matrix<3,X_DIM>& P, const Matrix<X_DIM,1>& xInit1, const Matrix<X_DIM,1>& xInit2, const int& radius) {
	Transform = !(P*G);
	Translate = -P*F*(xInit1-xInit2);
	for(int i = 0; i<NUM_POINTS; i++)
		ellipsoids.push_back(Transform*(points[i] + Translate));
}

// Finds the points of the ellipsoid obstacle within the reachable velocity space
inline void findReachableObstacle(const std::vector< Matrix<3,1> >& ellipsoids, std::vector< Matrix<3,1> >& reachablePoints, const State& x0, const State& x1, const double& maxSpeed) {
	// Assumes reachable velocities are an ellipsoid
	double xr, yr, zr;
	double xc, yc, zc;
	xr = maxSpeed;
	yr = maxSpeed;
	zr = maxSpeed;
	xc = x0[3]-x1[3];
	yc = x0[4]-x1[4];
	zc = x0[5]-x1[5];

	for(int i = 0; i<ellipsoids.size(); i++) {
		if( (pow((ellipsoids[i][0]-xc),2)/pow(xr,2) + pow((ellipsoids[i][1]-yc),2)/pow(yr,2) + pow((ellipsoids[i][2]-zc),2)/pow(zr,2)) < 1.0 )
			reachablePoints.push_back(ellipsoids[i]);
	}
}


// Should be replaced with the GJKEPA implementation in Bullet which also finds penetration depth if the piont is inside - 4/2/13
// Runs the GJK algorithm to find the distance between a given point (relative velocity) and a convex hull (LQRij intersection with reachable velocity)
// GJK only provides output NOT equal to 0 when the point is outside the hull, returns 0 if point is within hull
inline void run_gjk(const State& x0, const State& x1, const std::vector<Matrix<3,1>>& reachablePoints, double& distance, Matrix<3,1>& normal) {
	struct Object_structure VrelPoint;
	VrelPoint.numpoints = 1;
	
	REAL points1[1][3];
	points1[0][0] = x0[3]-x1[3];
	points1[0][1] = x0[4]-x1[4];
	points1[0][2] = x0[5]-x1[5];
	VrelPoint.vertices = points1;
	int ring1[3];
	ring1[0] = 1;
	ring1[1] = 0;
	ring1[2] = -1;
	VrelPoint.rings = ring1;
								
	struct Object_structure ConvexHull;
	ConvexHull.numpoints = reachablePoints.size();
	REAL points2[NUM_POINTS*OBSTACLE_STEPS][3];
	for(int m = 0; m < reachablePoints.size(); m++) {
		for(int n = 0; n < 3; n++) {
			points2[m][n] = reachablePoints[m][n];
		}
	}
	ConvexHull.vertices = points2;
	ConvexHull.rings = NULL;
	
	REAL wptHull[3];
	REAL wptVrel[3];
	wptHull[0] = wptHull[1] = wptHull[2] = wptVrel[0] = wptVrel[1] = wptVrel[2] = 0;
	distance = sqrt(gjk_distance(&VrelPoint, NULL, &ConvexHull, NULL, wptVrel, wptHull, NULL, 0)); // GJK returns distance^2
	
	normal[0] = ((wptVrel[0]-wptHull[0])/(distance));
	normal[1] = ((wptVrel[1]-wptHull[1])/(distance));
	normal[2] = ((wptVrel[2]-wptHull[2])/(distance));
}

// Returns bool to determine if a relative velocity is within the LQRij or not
void pointInHull(const std::vector<Matrix<3,1>>& reachablePoints, bool& insideHull, const State& x0, const State& x1) {
	double distance;
	Matrix<3,1> normal;
	double epsilon = 0.0001;
	run_gjk(x0,x1,reachablePoints, distance, normal);
	if(distance < epsilon && distance > -1*epsilon)
		insideHull = true;
	else
		insideHull = false;
}

// Needs to be eliminated using Bullet for convex hull - 4/2/13
// Writes necessary values to the text files for external running of the convexhull exe 
void convexHullFiles(const std::vector<Matrix<3,1>>& reachablePoints, std::vector<Matrix<4,1>>& planeEq, std::vector<Matrix<3,1>>& pointOnPlane) {
	std::ofstream outputFile;
	outputFile.open("pointList.txt");
	outputFile << "3" << std::endl << reachablePoints.size() << std::endl;
	for (int i = 0; i < reachablePoints.size(); i++)
		outputFile << reachablePoints[i][0] << " " << reachablePoints[i][1] << " " << reachablePoints[i][2] << std::endl;
	outputFile.close();
	system("qconvex n TO \"Planes.txt\" < \"pointList.txt\"");
	system("qconvex Fv TO \"facetVertices.txt\" < \"pointList.txt\"");
	std::string dump, first, second, third, fourth;	// Strings for reading values from text file
	int count = 0;						// dummy counter to skip first line of text read
	std::ifstream inputFile;
	inputFile.open("Planes.txt");
	int dim, numFacets;
	inputFile >> dim;
	inputFile >> numFacets;
	Matrix<4,1> tempPlane;				// Temp plane equation for pushback of planeEq
	for (int i = 0; i < numFacets; ++i) {
		inputFile >> tempPlane[0];
		inputFile >> tempPlane[1];
		inputFile >> tempPlane[2];
		inputFile >> tempPlane[3];
		planeEq.push_back(tempPlane);
	}
	inputFile.close();

	Matrix<3,1> tempPoint;
	int index;
	inputFile.open("facetVertices.txt");
	count = 0;
	inputFile >> numFacets;
	for (int i = 0; i < numFacets; ++i) {
		int numVertices;
		inputFile >> numVertices;
		for (int j = 0; j < numVertices; ++j) {
			int vertexId;
			inputFile >> vertexId;
			if (j == 0) {
				tempPoint[0] = reachablePoints[vertexId][0]; tempPoint[1] = reachablePoints[vertexId][1]; tempPoint[2] = reachablePoints[vertexId][2];
				pointOnPlane.push_back(tempPoint);
			}
		}
	}

}

// Find the convex hull of the reachable points
void convexHull(const std::vector<Matrix<3,1>>& reachablePoints, const State& x0, const State& x1, double& distance, Matrix<3,1>& normal) {
	// Write the reachable pioints to a file to input into qconvex.exe
	// Use GJK, and if distance is 0, then the point is inside the hull.  If the point is inside the hull, the iterative distance method
	//     must be used to find the shortest vector out of the obstacle.
	// If the point is outside the hull, GJK can be used to find the shortest distance to the hull.
	system("qconvex n TO \"dummyFile.txt\" < \"pointList.txt\"");
	system("qconvex Fv TO \"dummyFile.txt\" < \"pointList.txt\"");
	Matrix<3,1> relVel;
	relVel[0] = x0[3]-x1[3]; relVel[1] = x0[4]-x1[4]; relVel[2] = x0[5]-x1[5];  // Relative velocity term

	std::vector<Matrix<4,1>> planeEq;	// Vector of plane equations	
	std::vector<Matrix<3,1>> pointOnPlane;

	convexHullFiles(reachablePoints, planeEq, pointOnPlane);

	int dummy = 0;
	distance = abs(planeEq[0][0]*(relVel[0]-pointOnPlane[0][0]) + 
			   planeEq[0][1]*(relVel[1]-pointOnPlane[0][1]) + 
			   planeEq[0][2]*(relVel[2]-pointOnPlane[0][2]));
	double tempDist;
	for(int i = 1; i < pointOnPlane.size(); i++) {
		tempDist = abs(planeEq[i][0]*(relVel[0]-pointOnPlane[i][0]) + 
				   planeEq[i][1]*(relVel[1]-pointOnPlane[i][1]) + 
			       planeEq[i][2]*(relVel[2]-pointOnPlane[i][2]));
		if(tempDist < distance) {
			distance = tempDist;
			normal[0] = planeEq[i][0]; normal[1] = planeEq[i][1]; normal[2] = planeEq[i][2];
		}
	}
}

// Error value for ORCA library functions
static const float RVO_EPSILON = 0.00001f;

// Plane struct used for half-planes
struct Plane
{
	/*!
	 *  @brief		A point on the plane.
	 */
	Vector3 point;

	/*!
	 *  @brief		The normal to the plane.
	 */
	Vector3 normal;
};

// Line struct for 2-D half-planes
struct Line
{
	/*!
	 *  @brief		A point on the directed line.
	 */
	Vector3 point;

	/*!
	 *  @brief		The direction of the directed line.
	 */
	Vector3 direction;
};

// From ORCA_3D
bool linearProgram1(const std::vector<Plane>& planes, size_t planeNo, Line& line, float radius, const Vector3& optVelocity, bool directionOpt, Vector3& result)
  {
    const float dotProduct = line.point * line.direction;
    const float discriminant = sqr(dotProduct) + sqr(radius) - absSq(line.point);

    if (discriminant < 0.0f) {
      /* Max speed sphere fully invalidates line. */
      return false;
    }

    const float sqrtDiscriminant = std::sqrt(discriminant);
    float tLeft = -dotProduct - sqrtDiscriminant;
    float tRight = -dotProduct + sqrtDiscriminant;

    for (size_t i = 0; i < planeNo; ++i) {
      const float numerator = (planes[i].point - line.point) * planes[i].normal;
      const float denominator = line.direction * planes[i].normal;
      
      if (sqr(denominator) <= RVO_EPSILON) {
        /* Lines line is (almost) parallel to plane i. */
        if (numerator > 0.0f) {
          return false;
        } else {
          continue;
        }
      }

      const float t = numerator / denominator;

      if (denominator >= 0.0f) {
        /* Plane i bounds line on the left. */
        tLeft = max(tLeft, t);
		// tLeft = std::max(tLeft,t);
      } else {
        /* Plane i bounds line on the right. */
        tRight = min(tRight, t);
		//tRight = std::min(tRight,t);
      }

      if (tLeft > tRight) {
        return false;
      }
    }

    if (directionOpt) {
      /* Optimize direction. */
      if (optVelocity * line.direction > 0.0f) {
        /* Take right extreme. */
        result = line.point + tRight * line.direction;
      } else {
        /* Take left extreme. */
        result = line.point + tLeft * line.direction;
      }
    } else {
      /* Optimize closest point. */
      const float t = line.direction * (optVelocity - line.point);

      if (t < tLeft) {
        result = line.point + tLeft * line.direction;
      } else if (t > tRight) {
        result = line.point + tRight * line.direction;
      } else {
        result = line.point + t * line.direction;
      }
    }

    return true;
  }

// From ORCA_3D
  bool linearProgram2(const std::vector<Plane>& planes, size_t planeNo, float radius, const Vector3& optVelocity, bool directionOpt, Vector3& result)
  {
    const float planeDist = planes[planeNo].point * planes[planeNo].normal;
    const float planeDistSq = sqr(planeDist);
    const float radiusSq = sqr(radius);
    if (planeDistSq > radiusSq) {
      /* Max speed sphere fully invalidates plane planeNo. */
      return false;
    }
    
    const float planeRadiusSq = radiusSq - planeDistSq;
    const Vector3 planeCenter = planeDist * planes[planeNo].normal;

    if (directionOpt) {
      /* Project direction optVelocity on plane planeNo */
      const Vector3 planeOptVelocity = optVelocity - (optVelocity * planes[planeNo].normal) * planes[planeNo].normal;
      const float planeOptVelocityLengthSq = absSq(planeOptVelocity); 
      if (planeOptVelocityLengthSq <= RVO_EPSILON) {
        result = planeCenter;
      } else {
        result = planeCenter + std::sqrt(planeRadiusSq / planeOptVelocityLengthSq) * planeOptVelocity; 
      }
    } else {
      /* Project point optVelocity on plane planeNo */
      result = optVelocity + ((planes[planeNo].point - optVelocity) * planes[planeNo].normal) * planes[planeNo].normal;
      /* If outside planeCircle, project on planeCircle */
      if (absSq(result) > radiusSq) {
        const Vector3 planeResult = result - planeCenter;
        const float planeResultLengthSq = absSq(planeResult);
        result = planeCenter + std::sqrt(planeRadiusSq / planeResultLengthSq) * planeResult;
      }
    }

    for (size_t i = 0; i < planeNo; ++i) {
      if (planes[i].normal * (planes[i].point - result) > 0.0f) {
        /* Result does not satisfy constraint i. Compute new optimal result. */

        /* Compute intersection line of plane i and plane planeNo. */
        Vector3 crossProduct = cross(planes[i].normal, planes[planeNo].normal);
      
        if (absSq(crossProduct) <= RVO_EPSILON) {
          /* Planes planeNo and i are (almost) parallel, and
           * plane i fully invalidates plane planeNo. 
           */
          return false;
        }

        Line line;
        line.direction = normalize(crossProduct);
        const Vector3 lineNormal = cross(line.direction, planes[planeNo].normal);
        line.point = planes[planeNo].point + (((planes[i].point - planes[planeNo].point) * planes[i].normal) / (lineNormal * planes[i].normal)) * lineNormal;

        if (!linearProgram1(planes, i, line, radius, optVelocity, directionOpt, result)) {
          return false;
        }
      }
    }
    
    return true;
  }

  // From ORCA_3D
  size_t linearProgram3(const std::vector<Plane>& planes, double radius, const Vector3& optVelocity, bool directionOpt, Vector3& result)
  {
    if (directionOpt) {
      /* 
       * Optimize direction. Note that the optimization velocity is of unit
       * length_ in this case.
       */
      result = optVelocity * radius;
    } else if (absSq(optVelocity) > sqr(radius)) {
      /* Optimize closest point and outside circle. */
      result = normalize(optVelocity) * radius;
    } else {
      /* Optimize closest point and inside circle. */
      result = optVelocity;
    }

    for (size_t i = 0; i < planes.size(); ++i) {
      if (planes[i].normal * (planes[i].point - result) > 0.0f) {
        /* Result does not satisfy constraint i. Compute new optimal result. */
        const Vector3 tempResult = result;
        
        if (!linearProgram2(planes, i, radius, optVelocity, directionOpt, result)) {
          result = tempResult;
          return i;
        }
      }
    }

    return planes.size();
  }

  // From ORCA_3D
  void linearProgram4(const std::vector<Plane>& planes, size_t beginPlane, float radius, Vector3& result)
  {
    float distance = 0.0f;

    for (size_t i = beginPlane; i < planes.size(); ++i) {
      if (planes[i].normal * (planes[i].point - result) > distance) {
        // Result does not satisfy constraint of plane i.
        std::vector<Plane> projPlanes;

        for (size_t j = 0; j < i; ++j) {
          Plane plane;

          const Vector3 crossProduct = cross(planes[j].normal, planes[i].normal);

          if (absSq(crossProduct) <= RVO_EPSILON) {
            // Plane i and plane j are (almost) parallel.
            if (planes[i].normal * planes[j].normal > 0.0f) {
              // Plane i and plane j point in the same direction.
              continue;
            } else {
              // Plane i and plane j point in opposite direction.
              plane.point = 0.5f * (planes[i].point + planes[j].point);
            }
          } else {
            // Plane.point is point on line of intersection between plane i and plane j
            const Vector3 lineNormal = cross(crossProduct, planes[i].normal);
            plane.point = planes[i].point + (((planes[j].point - planes[i].point) * planes[j].normal) / (lineNormal * planes[j].normal)) * lineNormal;
          }

          plane.normal = normalize(planes[j].normal - planes[i].normal);
          projPlanes.push_back(plane);
        }

        const Vector3 tempResult = result;

        if (linearProgram3(projPlanes, radius, planes[i].normal, true, result) < projPlanes.size()) {
          // This should in principle not happen.  The result is by definition
          // already in the feasible region of this linear program. If it fails,
          // it is due to small floating point error, and the current result is
          // kept.
          result = tempResult;
        }

        distance = planes[i].normal * (planes[i].point - result);
      }
    }
  }

  // Create a halfplane RCA_{ij} given a relative velocity (in or out of convex hull) and it's distance and normal to intersect (point out of hull) or escape (point in hull)
  void createHalfPlanes(const State& x, double& distance, Matrix<3,1>& normal, const bool& insideHull, std::vector<Plane>& orcaPlanes_) {
	Plane tempPlane;
	Vector3 tempNormal(normal[0], normal[1], normal[2]);
	double mult;
	if(insideHull)
		mult = 1.0;
	else
		mult = -1.0;

	Vector3 tempVector(x[3]+mult*distance*normal[0], x[4]+mult*distance*normal[1], x[5]+mult*distance*normal[2]); 
	tempPlane.normal = tempNormal;
	tempPlane.point = tempVector;
	orcaPlanes_.push_back(tempPlane);
  }

  // Find a new velocity that is as close as possible to goal velocity using linear programming on a set of half-planes
  void calculateNewV(std::vector<Plane>& orcaPlanes_, const Matrix<3,1>& goalVel, Matrix<3,1>& newVel) {
	double maxSpeed_ = 100;
	Vector3 prefVelocity_(goalVel[0], goalVel[1], goalVel[2]);
	Vector3 newVelocity_;

	size_t planeFail = linearProgram3(orcaPlanes_, maxSpeed_, prefVelocity_, false, newVelocity_);
    if (planeFail < orcaPlanes_.size())
      linearProgram4(orcaPlanes_, planeFail, maxSpeed_, newVelocity_);

	newVel[0] = newVelocity_.x(); newVel[1] = newVelocity_.y(); newVel[2] = newVelocity_.z();
	orcaPlanes_.clear();
  }

  // Main simulation func
int _tmain(int argc, _TCHAR* argv[]) {
	// Initialize state setup
	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,U_DIM> B;
	Matrix<X_DIM,1> c;

	// Set up simulation/visualization
	setup();
	std::vector< Quadrotor* > qlist(NUM_QUADS);

	// Define goal thrust
	Input uGoal;
	uGoal[0] = uGoal[1] = uGoal[2] = uGoal[3] = nominalInput;

	// LQR path penalty matrices
	Qx = zeros<X_DIM,X_DIM>();
	Qv = zeros<V_DIM,V_DIM>();
	Qv(0,0) = Qv(1,1) = Qv(2,2) = 100;
	Qp = zeros<V_DIM,V_DIM>();
	Qp(0,0) = Qp(1,1) = Qp(2,2) = 0.1;

	// LQR input penalty
	R = 5*identity<U_DIM>();

	// Motion noise variances
	M = 0.000000001 * identity<X_DIM>(); // motion noise variance
	N = 0.000000001 * identity<Z_DIM>(); // observation noise variance

	// Begin setup for quadrotors
	for( int i = 0; i < NUM_QUADS; i++) {
		qlist[i] = new Quadrotor();
		qlist[i]->setupQuadVisualization();
	}

	// Initalize a quad state
	State xInit = zeros<X_DIM>();
	xInit[12] = xInit[13] = xInit[14] = xInit[15] = nominalInput;

	Rotation RotInit = identity<3>();
	Matrix<X_DIM,X_DIM> Pinit = 0.000000001 * identity<X_DIM>();

	State xGoal = zeros<X_DIM>();
	xGoal[12] = xGoal[13] = xGoal[14] = xGoal[15] = nominalInput;

	std::vector<Position> pGoal;
	Position tempPGoal;
	Rotation RotGoal = identity  <3>();
	
	// Initialize quadrotor beginning and goal positions
	double height=0.0;
	double diag = std::sqrt(3.0);
	Matrix<3,1> One, Two, Three, Four, Five, Six;
	One[0] = 6.0; One[1] = 0.0; One[2] = height;
	Two[0] = 3.0; Two[1] = 3.0*diag; Two[2] = height;
	Three[0] = -3.0; Three[1] = 3.0*diag; Three[2] = height;
	Four[0] = -6.0; Four[1] = 0.0; Four[2] = height;
	Five[0] = -3.0; Five[1] = -3.0*diag; Five[2] = height;
	Six[0] = 3.0; Six[1] = -3.0*diag; Six[2] = height;

	// Quadrotor 1
	xInit[0] = One[0]; xInit[1] = One[1]; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0.0;
	tempPGoal[0] = Four[0]; tempPGoal[1] = Four[1]; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[0]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[0], uGoal, RotGoal);
	if(NUM_QUADS > 1) {
		// Quadrotor 2
		xInit[0] = Two[0]; xInit[1] = Two[1]; xInit[2] = height;
		xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0.0;
		tempPGoal[0] = Five[0]; tempPGoal[1] = Five[1]; tempPGoal[2] = height;
		pGoal.push_back(tempPGoal);
		qlist[1]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[1], uGoal, RotGoal); 
	}
	if(NUM_QUADS > 2) {
		// Quadrotor 3
		xInit[0] = Three[0]; xInit[1] = Three[1]; xInit[2] = height;
		xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0.0;
		tempPGoal[0] = Six[0]; tempPGoal[1] = Six[1]; tempPGoal[2] = height;
		pGoal.push_back(tempPGoal);
		qlist[2]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[2], uGoal, RotGoal);
	}
	if(NUM_QUADS > 3) {
		// Quadrotor 4
		xInit[0] = Four[0]; xInit[1] = Four[1]; xInit[2] = height;
		xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0.0;
		tempPGoal[0] = One[0]; tempPGoal[1] = One[1]; tempPGoal[2] = height;
		pGoal.push_back(tempPGoal);
		qlist[3]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[3], uGoal, RotGoal);
	}
	if(NUM_QUADS > 4) {
		// Quadrotor 5
		xInit[0] = Five[0]; xInit[1] = Five[1]; xInit[2] = height;
		xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
		tempPGoal[0] = Two[0]; tempPGoal[1] = Two[1]; tempPGoal[2] = height;
		pGoal.push_back(tempPGoal);
		qlist[4]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[4], uGoal, RotGoal);
	}
	if(NUM_QUADS > 5) {
		// Quadrotor 6
		xInit[0] = Six[0]; xInit[1] = Six[1]; xInit[2] = height;
		xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
		tempPGoal[0] = Three[0]; tempPGoal[1] = Three[1]; tempPGoal[2] = height;
		pGoal.push_back(tempPGoal);
		qlist[5]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[5], uGoal, RotGoal);
	}
	/*// Quadrotor 6
	xInit[0] = 0.0; xInit[1] = 4.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = 0.0; tempPGoal[1] = -4.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[6]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[6], uGoal, RotGoal);
	//Quadrotor 7
	xInit[0] = 0.0; xInit[1] = -4.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = 0.0; tempPGoal[1] = 4.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[7]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[7], uGoal, RotGoal);
	double arm  = sqrt(8.0);
	double arm2 = sqrt(18.0);
	double arm3 = sqrt(2.0);
	//Quadrotor 8
	xInit[0] = arm; xInit[1] = arm; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm2; tempPGoal[1] = -arm2; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[8]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[8], uGoal, RotGoal);
	//Quadrotor 9
	xInit[0] = arm; xInit[1] = -arm; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm2; tempPGoal[1] = arm2; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[9]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[9], uGoal, RotGoal);
	//Quadrotor 10
	xInit[0] = -arm; xInit[1] = -arm; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm2; tempPGoal[1] = arm2; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[10]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[10], uGoal, RotGoal);
	//Quadrotor 11
	xInit[0] = -arm; xInit[1] = arm; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm2; tempPGoal[1] = -arm2; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[11]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[11], uGoal, RotGoal);
	//Quadrotor 12
	xInit[0] = arm2; xInit[1] = arm2; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm; tempPGoal[1] = -arm; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[12]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[12], uGoal, RotGoal);
	//Quadrotor 13
	xInit[0] = arm2; xInit[1] = -arm2; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm; tempPGoal[1] = arm; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[13]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[13], uGoal, RotGoal);
	//Quadrotor 14
	xInit[0] = -arm2; xInit[1] = -arm2; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm; tempPGoal[1] = arm; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[14]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[14], uGoal, RotGoal);
	//Quadrotor 15
	xInit[0] = -arm2; xInit[1] = arm2; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm; tempPGoal[1] = -arm; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[15]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[15], uGoal, RotGoal);
	//Quadrotor 16
	xInit[0] = 2.0; xInit[1] = 0.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -8.0; tempPGoal[1] = 0.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[16]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[16], uGoal, RotGoal);
	//Quadrotor 17
	xInit[0] = -2.0; xInit[1] = 0.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = 8.0; tempPGoal[1] = 0.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[17]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[17], uGoal, RotGoal);
	//Quadrotor 18
	xInit[0] = 0.0; xInit[1] = 2.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = 0.0; tempPGoal[1] = -8.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[18]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[18], uGoal, RotGoal);
	//Quadrotor 19
	xInit[0] = 0.0; xInit[1] = -2.0; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = 0.0; tempPGoal[1] = 8.0; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[19]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[19], uGoal, RotGoal);
		//Quadrotor 20
	xInit[0] = arm3; xInit[1] = arm3; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm3; tempPGoal[1] = -arm3; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[20]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[20], uGoal, RotGoal);
	//Quadrotor 21
	xInit[0] = arm3; xInit[1] = -arm3; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = -arm3; tempPGoal[1] = arm3; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[21]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[21], uGoal, RotGoal);
	//Quadrotor 22
	xInit[0] = -arm3; xInit[1] = -arm3; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm3; tempPGoal[1] = arm3; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[22]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[22], uGoal, RotGoal);
	//Quadrotor 23
	xInit[0] = -arm3; xInit[1] = arm3; xInit[2] = height;
	xInit[3] = 0.0; xInit[4] = 0.0; xInit[5] = 0;
	tempPGoal[0] = arm3; tempPGoal[1] = -arm3; tempPGoal[2] = height;
	pGoal.push_back(tempPGoal);
	qlist[23]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, pGoal[23], uGoal, RotGoal);*/
	
	/*for(int i = 0; i < NUM_QUADS; ++i) {
		xInit[0] = -5+random()*10; xInit[1] = -5+random()*10; xInit[2] = -5+random()*10;
		tempPGoal[0] = -5+random()*10; tempPGoal[1] = -5+random()*10; tempPGoal[2] = -5+random()*10;
		pGoal.push_back(tempPGoal);
		qlist[i]->setupQuadrotors(xInit,RotInit,Pinit, xGoal, pGoal[i], uGoal, RotGoal);
	}
	*/

	// Configuration matrix mapping position from state
	Matrix<3,X_DIM> C = zeros<3,X_DIM>();
	C(0,0) = C(1,1) = C(2,2) = 1;
	Matrix<V_DIM,X_DIM> V = zeros<V_DIM,X_DIM>();
	V(0,3) = V(1,4) = V(2,5) = 1;

	// Setup discretized FG matrices
	Matrix<X_DIM,X_DIM> Ft = identity<X_DIM>();		// Ft = A^t
	Matrix<X_DIM,V_DIM> Gt = zeros<X_DIM,V_DIM>();  // Gt = (A^(t-1) + A^(t-2) + A^(t-3) ... + A + 1)B

	CAL_DestroyGroup(rotordummy);				// erase original quad visualizer that was copied
	CAL_DestroyGroup(ellipsedummy);				// erase bounding ellipse that was copied

	for (size_t i = 0; i < NUM_QUADS; ++i) {	// Loop through all quads
		qlist[i]->findMatrices(A,B,c);			// Linearize/discretize system
		qlist[i]->vGoal = qlist[i]->findVGoal(); // Set goal velocities
	}

	std::vector< Matrix<3,1> > points;			// Vector to save the points of a sampled sphere
	createSpheres(points);						// Function that samples the sphere
	bool insideHull = false;					// Bool to test if the relative velocity point is inside or outside the convex hull
	int numVertices;							// Number of vertices of the convex hull
	Matrix<4,1> halfPlaneVector;				// Elements 0-2 = normal vector out of the convex hull, element 3 = magnitude of distance vector
	double distance;							// Distance for halfPlaneVector;
	std::vector<double> distances;				// Vector of distances to use in linear programming
	Matrix<3,1> normalVector;					// Normal vector for halfPlaneVector
	std::vector<Matrix<3,1>> normalVectors;		// Vector of normal vectors, also for linear programming
	double maxSpeed = 30.0;						// Maximum velocity quadrotor can obtain
	int numObstacle;							// Number of points in the intersection of reachable velocity set and LQR set
	std::vector<Plane> orcaPlanes_;				// Empty vector of halfplanes for linearprogramming
	std::cout << "Start!" << std::endl;			// Indicator that simulation has started
	for(size_t t = 0; t < TIME_STEPS; ++t) {	// Loop until simulation time has been reached
		std::cout << t << std::endl;			// Output simulation time while calculating
		for(size_t i = 0; i  < NUM_QUADS; ++i) { // Loop over all quadrotors
			distances.clear();					// clear the half-plane distances
			normalVectors.clear();				// clear the half-plane normals
			for(size_t j = 0; j < NUM_QUADS; ++j) { // Loop over all OTHER quads for a given quad
				qlist[i]->ellipsoids.clear();	// clear the ellipsoids for a given quad
				qlist[i]->reachablePoints.clear();	// clear reachable points
				if(i==j);						// skip comparing quad to itself
				else {
					Ft = identity<X_DIM>();		// initialize F[t]
					Gt = zeros<X_DIM,V_DIM>();	// initialize G[t]
					for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {	// loop over time tau
						findFG(A,B,qlist[i]->L, qlist[i]->E, Ft, Gt);	// find F[t] and G[t] for quad pair
						createObstacle(qlist[i]->ellipsoids, points, qlist[i]->Transform, qlist[i]->Translate, t, Ft, Gt, C, qlist[i]->x, qlist[j]->x,1);
						// ^^ Create the LQR obstacle for a given F[t], G[t], and relative state.
						// Only makes one ellipsoid at a time, but through the loop over tau creates the union of all ellipses,
						// thus building LQRij
					}
					findReachableObstacle(qlist[i]->ellipsoids, qlist[i]->reachablePoints, qlist[i]->x, qlist[j]->x, maxSpeed);
					// ^^ find points that are reacahble and within the LQRij
					numObstacle= qlist[i]->reachablePoints.size(); // Number of points in the obstacle. GJK errors with less than 4
					if(numObstacle > 4) {	
					    pointInHull(qlist[i]->reachablePoints, insideHull, qlist[i]->x, qlist[j]->x); // Sets bool for if point is within LQRij
					    if(insideHull){		// if vrel within LQRij
						    convexHull(qlist[i]->reachablePoints, qlist[i]->x, qlist[j]->x, distance, normalVector); // Run convex hull and get distance to nearest point
						}
					    else {				// if vrel outside LQRij
						    run_gjk(qlist[i]->x, qlist[j]->x, qlist[i]->reachablePoints, distance, normalVector); // run GJK to get distance to intersection
						}
					    distance *=0.5; // half the distance for the reciprocity aspect of algorithm
					    createHalfPlanes(qlist[i]->x, distance, normalVector, insideHull, orcaPlanes_); // Create a half-plane for a given vrel and LQRij
					}
				}
			}
			calculateNewV(orcaPlanes_, qlist[i]->vGoal, qlist[i]->newV); // Once all half-planes for one quad to all others is populated, find the new collision-avoiding vel.
		}
		for(size_t i = 0; i < NUM_QUADS; ++i) { // Once all goalvels are updated, propagate them through the control to the next time-step of simulation
			qlist[i]->vGoal = qlist[i]->newV;	// Update vGoal
			Input u = qlist[i]->findU();		// find LQR goal output
			qlist[i]->propagateU(u);			// send output to robot
			kalmanFilter1(qlist[i]->x, qlist[i]->Rot, u, M, qlist[i]->P);	// get state estimate with kalman
			Observation z = sampleGaussian(h(qlist[i]->xTrue, qlist[i]->RotTrue), N);
			kalmanFilter2(qlist[i]->x, qlist[i]->Rot, z, N, qlist[i]->P);
			qlist[i]->visualize(t*dt);			// update visualization
			qlist[i]->vGoal = qlist[i]->findVGoal();	// reset vGoal 
		}
	}
	std::cout << "Done" << std::endl; // Output message notifying simulation is done calculating and ready to view

	int k; 
	std::cin >> k; // Pause so user can use the callisto UI

	CAL_End(); // end callisto
	return 0;  // exit main
}
