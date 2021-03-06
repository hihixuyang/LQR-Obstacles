#include "simulator2.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "gjk.h"
#include "Vector3.h"
//#include <math.h>

#define NUM_QUADS 2			// Number of quadrotors
#define TIME_STEPS 300		// Time steps for total simulation at 1/30 seconds per step
#define OBSTACLE_STEPS 45	// Time to check for obstacle
#define NUM_POINTS 100		// Number of points to make up the ellipse
#define XYRADIUS 0.26
#define ZRADIUS 0.26 
bool drawEllipses = false;

clock_t start = 0, stop = 0;
clock_t fstart = 0, fstop = 0;
double duration=0, ftotal=0;

/*#define X_DIM 16			// State dimension
#define V_DIM 3				// Velocity control dimension
#define U_DIM 4				// Control input dimension
#define Z_DIM 6				// Observation dimension

typedef Matrix<X_DIM> State;
typedef Matrix<U_DIM> Input;
typedef Matrix<Z_DIM> Observation;
typedef Matrix<3,3>   Rotation;
typedef Matrix<3,1>   Velocity;*/

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
//std::vector< int > CALellipses((TIME_STEPS+2));
int rotordummy;
int ellipsedummy;
int cal_scene;
int cal_floor;
int cal_ellipse;
//std::vector< int > ellipses(NUM_POINTS);

int counter=0;

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
  //State xInit;
  State xGoal;
  Velocity vGoal;
  Position pGoal;
  Velocity newV;
  Input uGoal;
  Rotation RotGoal;
  
  int cal_quadrotor;
  int cal_ellipse;

  State x;
  Rotation Rot;
  State xTrue;
  Rotation RotTrue;
  Matrix<U_DIM,X_DIM> L;
  Matrix<U_DIM,V_DIM> E;
  Matrix<U_DIM,1> l; // ell*/
  Matrix<V_DIM,X_DIM> Lh;
  Matrix<V_DIM,V_DIM> Eh;

  Matrix<X_DIM,X_DIM> P;

  std::vector< Matrix<3,1> > ellipsoids;		// Vector of points making up the whole obstacle
  std::vector< Matrix<3,1> > reachablePoints;	// Vector of points from the obstacle that are reachable
  Matrix<3,3> Transform;						// Transformation matrix PG^-1
  Matrix<3,1> Translate;						// Translation Matrix
  std::vector< Matrix<4,1> > planes;			// Convex hull facets' plane equation coefficients
  std::vector< Matrix<3,1> > normals;			// Convex hull normal vector
  std::vector< Matrix<3,3> > vertex;			// Convex hull facet's plane vertex
  Matrix<4,1> shortestDistance;					// Escape vector. Elements 0-2 = direction, 3 = magnitude 

  void setupQuadrotors(const State& xinit, const Rotation& rotinit, const Matrix<X_DIM,X_DIM> Pinit, const State& xfinal, const Position& pfinal, const Input& uHover, const Rotation& rotfinal) {
    x = xinit;
    Rot = rotinit;
    xGoal = xfinal;
    pGoal = pfinal;
    uGoal = uHover;
    RotGoal = rotfinal;
    P = Pinit;

    xTrue = sampleGaussian(x,P);
    RotTrue = Rot;
    RotTrue = RotTrue * exp(skewSymmetric((xTrue).subMatrix<3,1>(6,0)));
    xTrue[6] = 0; xTrue[7] = 0; xTrue[8] = 0;
    //P = Pinit;
  }

  void setupQuadVisualization() {
    CAL_CloneGroup(&cal_quadrotor, rotordummy, 0, false, "Quadrotor");
    CAL_CloneGroup(&cal_ellipse, ellipsedummy, 0, false, "ellipse");
  }

  void visualize(double t) {
    float p[3] = {(float) xTrue[0], (float) xTrue[1], (float) xTrue[2]};
    Matrix<4> q = quatFromRot(RotTrue);
    float o[4] = {(float) q[0], (float) q[1], (float) q[2], (float) q[3]};
    CAL_AddGroupKeyState(cal_quadrotor, (float) t, p, o);
	CAL_AddGroupKeyState(cal_ellipse, (float) t, p);
    //CAL_SetGroupPosition(cal_quadrotor, p[0], p[1], p[2]);
    //CAL_SetGroupQuaternion(cal_quadrotor, q[0], q[1], q[2], q[3]);


    Matrix<3,3> V, E;
    jacobi(P.subMatrix<3,3>(0,0), V, E);
    q = quatFromRot(V);
    float pos[3] = {(float) x[0], (float) x[1], (float) x[2]};
    float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
    float scale[3] = {(float) (2*sqrt(E(0,0))), (float) (2*sqrt(E(1,1))), (float) (2*sqrt(E(2,2)))};
    //CAL_AddGroupKeyState(cal_ellipse, (float) t, pos, quat, scale);
	//CAL_SetGroupPosition(cal_ellipse, p[0], p[1], p[2]);
    //CAL_SetGroupQuaternion(cal_ellipse, q[0], q[1], q[2], q[3]);
  }

  void findMatrices(Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,1>& c) {
    controlMatrices(uGoal, xGoal, A, B, c, L, E, l, Lh, Eh);
  }

  Input findU() {
    return riccatiControllerSteady(x, Rot, vGoal, RotGoal, uGoal, L, E, l);
  }

   Velocity findVGoal() {
    return riccatiControllerSteadyPosition(x,Rot, pGoal, RotGoal, vGoal, Lh, Eh);
  }

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
  length      = 0.3429/2;          // distance between center and rotor, m
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

  double volumeWidth   = 8.0; // m
  double volumeLength  = 8.0;  // m
  double volumeHeight  = 3.3528;  // m

// visualization
  CAL_Initialisation(true, true, true);
	int obj;
  
  // Quadrotor
  CAL_CreateGroup(&rotordummy, 0, false, "QuadRotor");
  CAL_SetGroupColor(rotordummy, 0.05, 0.05, 0.05);
  CAL_CreateBox(rotordummy, 2*length, beamWidth, beamHeight, 0, 0, 0);
  CAL_CreateBox(rotordummy, beamWidth, 2*length, beamHeight, 0, 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, length, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, -length, 0, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, 0, length, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, motorRadius, motorHeight, 0, -length, beamHeight / 2 + motorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, length, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, -length, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, 0, length, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, beamRadius, beamHeight, 0, -length, 0, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, -length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, 0, length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateCylinder(rotordummy, rotorRadius, rotorHeight, 0, -length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
  CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
  CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
  CAL_CreateBox(rotordummy, centerSide, centerSide, beamHeight, 0, 0, 0, &obj);
  CAL_SetObjectOrientation(obj, 0, 0, (float) (M_PI*0.25));
  CAL_CreateBox(rotordummy, flagLength, beamWidth + 0.001, beamHeight + 0.001, length / 1.65, 0, 0, &obj);
  CAL_SetObjectColor(obj, 1, 0.15, 0);

  float flagTriangle[18] = {length / 1.65 - flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
                            length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
                            length / 1.65 - flagLength / 2, 0, -beamHeight / 2};
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
  CAL_CreateGroup(&cal_floor, 0, false, "Floor");
  CAL_LoadTexture(0, "floor.png", 1);
  CAL_CreateBox(cal_floor, volumeWidth, volumeLength, 0.1, 0, 0, -0.05, &obj);
  CAL_SetObjectTexture (obj, 0, volumeWidth/tileSize, volumeLength/tileSize);

  // Ellipse
  CAL_CreateGroup (&ellipsedummy, 0, false, "Ellipse");
  CAL_SetGroupColor(ellipsedummy, 0, 1, 0, 0.5);
  CAL_CreateSphere(ellipsedummy, 1, 0, 0, 0, &obj);
  CAL_SetObjectScaling(obj, XYRADIUS, XYRADIUS, ZRADIUS);
  
  // Volume
  CAL_CreateGroup(&cal_scene, 0, false, "Scene");
  }

// Function to create a sphere(ellipse) in Callisto 4
/*void createEllipse(int& ellipse) {
	CAL_CreateGroup(&ellipse,0, false, "ellipse");
	//CAL_CreateSphere(ellipse, 1, 0,0,0);
	CAL_CreateSphere(ellipse, 1, 0,0,0);
}*/

// Performs a transformation on the ellipse by some matrix M
/*void transformEllipse(int& ellipse, const Matrix<3,3>& M) {
	CAL_matrix3 mat;
	// Changed to transpose of the matrix to check... May or may not be right.  Need to check
	mat[0][0] = M(0,0); mat[0][1] = M(0,1); mat[0][2] = M(0,2);
	mat[1][0] = M(1,0); mat[1][1] = M(1,1); mat[1][2] = M(1,2);
	mat[2][0] = M(2,0); mat[2][1] = M(2,1); mat[2][2] = M(2,2);
	CAL_SetGroupMatrix(ellipse, mat);
}*/

// Translates the elllipse to some position Pos
/*void translateEllipse(int& ellipse, const Matrix<3,1>& Pos) {
	CAL_scalar x; CAL_scalar y; CAL_scalar z;
	x = Pos[0]; y = Pos[1]; z = Pos[2];
	CAL_SetGroupPosition(ellipse,x,y,z);
}*/

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

/*// Generate random vector from normal distribution N(mean, var)
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
}*/

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

// Added V, changed Q to Qv, added Qx, and x^*
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
  // Matrix<X_DIM,X_DIM> S = Q;
  Matrix<X_DIM, X_DIM> S = ~V*Qv*V; // + Qx;
  //Matrix<X_DIM, V_DIM> T = ~V*Qv;
  Matrix<X_DIM, V_DIM> T = ~V*Qv;
  for (int i = 0; i < 300; ++i) {
     T = ~V*Qv + ~A*T - ~A*S*B*!(R + ~B*S*B)*~B*T;
     S = ~V*Qv*V + Qx + ~A*S*A - ~A*S*B*!(R + ~B*S*B)*(~B*S*A);
  }

   //Matrix<X_DIM,1>   v = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * (~A*S*B*!(R + ~B*S*B)*~B*S - ~A*S)*c;
   Matrix<X_DIM,1>     a = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>())*(Qx*xstar - ~A*S*c + ~A*S*B*!(R+~B*S*B)*~B*S*c);

   // feedback law: u = L*xTilde + E*xTildeGoal + l
   L = -!(R + ~B*S*B)*~B*S*A;
   E = !(R + ~B*S*B)*~B*T;
   l = -!(R + ~B*S*B)*(~B*S*c + ~B*a);

   Matrix<X_DIM,X_DIM> Qptilde = ~Pmap*Qp*Pmap + ~L*Rp*L;
   Matrix<V_DIM,V_DIM> Rtilde = ~E*Rp*E;
   Matrix<V_DIM,X_DIM> Ptilde = ~E*Rp*L;
   Matrix<X_DIM,X_DIM> Atilde = A+B*L;
   Matrix<X_DIM,V_DIM> Btilde = B*E;

     // Solve riccati equation to find S (naive implementation)
  // Matrix<X_DIM,X_DIM> S = Q;
  Matrix<X_DIM, X_DIM> Stilde = ~Pmap*Qp*Pmap; // + Qx;
  Matrix<X_DIM, V_DIM> Ttilde = ~Pmap*Qp;
  for (int i = 0; i < 300; ++i) {
     Ttilde = ~Pmap*Qp + ~Atilde*Ttilde - (~Ptilde + ~Atilde*Stilde*Btilde)*!(Rtilde + ~Btilde*Stilde*Btilde)*~Btilde*Ttilde;
     Stilde = Qptilde + ~Atilde*Stilde*Atilde - (~Ptilde + ~Atilde*Stilde*Btilde)*!(Rtilde+~Btilde*Stilde*Btilde)*(Ptilde+~Btilde*Stilde*Atilde);
  }
//Matrix<X_DIM,1>   v = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * (~A*S*B*!(R + ~B*S*B)*~B*S - ~A*S)*c;
   //Matrix<X_DIM,1>     atilde = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>())*(Qx*xstar - ~A*S*c + ~A*S*B*!(R+~B*S*B)*~B*S*c);
  Matrix<X_DIM,1> atilde = zeros<X_DIM,1>();

   // feedback law: u = L*xTilde + E*xTildeGoal + l
   //Lh = -(!(Rtilde + ~Btilde*Stilde*Btilde)*(~Ptilde + (~Btilde*Stilde*Atilde)));
   Lh = -(!(Rtilde+~Btilde*Stilde*Btilde))*(Ptilde+~Btilde*Stilde*Atilde);
   Eh = -(!(Rtilde + ~Btilde*Stilde*Btilde))*(~Btilde*Ttilde);
}

// Needs to be updated to control off of vstar not xstar
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

inline Velocity riccatiControllerSteadyPosition(const State& x, const Rotation& R0, const Position& pGoal, const Rotation& RGoal, const Velocity& vGoal, const Matrix<V_DIM,X_DIM>& Lh, const Matrix<V_DIM,V_DIM>& Eh) {
  
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
  xTilde.insert(3,0, ~RLocal * (x.subMatrix<3,1>(3,0) - vGoal));
  xTilde.insert(6,0, errFromRot(~RLocal*R0)); //xTilde.insert(6,0, -(~R0*axis));
  xTilde.insert(9,0, x.subMatrix<3,1>(9,0));
  xTilde.insert(12,0, x.subMatrix<4,1>(12,0));
  
  return vGoal + Lh*xTilde + Eh*pTildeGoal;
}
/*inline Input riccatiControllerCurrent(const State& x, const Rotation& R0, const State& xGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<X_DIM, X_DIM>& Q, const Matrix<U_DIM, U_DIM>& R, Matrix<X_DIM, X_DIM>& S) {
	/*Matrix<X_DIM,X_DIM> A;
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
}*/

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


inline void findFG(const Matrix<X_DIM,X_DIM>& A, const Matrix<X_DIM,U_DIM>& B, const Matrix<U_DIM,X_DIM>& L, const Matrix<U_DIM,V_DIM>& E, Matrix<X_DIM,X_DIM>& F, Matrix<X_DIM,V_DIM>& G) {
	Matrix<X_DIM,X_DIM> Atilde;
	Atilde = A+(B*L);

	Matrix<X_DIM,V_DIM> Btilde;
	Btilde = (B*E);

	F = Atilde*F;
	G = Atilde*G + Btilde;
}

inline void findDistance(const State& x0, const State& x1, const std::vector<Matrix<3,1>>& reachablePoints, double& distance, Matrix<3,1>& normal){
	
}

// Creates a sphere of N points
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

// Takes the sphere of N points and makes an ellipse out of it.
inline void createObstacle(std::vector< Matrix<3,1> >& ellipsoids, const std::vector< Matrix<3,1> >& points, Matrix<3,3>& Transform, Matrix<3,1>& Translate, const size_t& t, const Matrix<X_DIM,X_DIM>& F, const Matrix<X_DIM,V_DIM>& G, const Matrix<3,X_DIM>& P, const Matrix<X_DIM,1>& xInit1, const Matrix<X_DIM,1>& xInit2) {
	Transform = !(P*G);
	//std::cout << Transform << std::endl << std::endl;
	Translate = -P*F*(xInit1-xInit2);
	//std::cout << Translate << std::endl << std::endl;
	for(int i = 0; i<NUM_POINTS; i++)
		ellipsoids.push_back(Transform*(points[i] + Translate));
	
	/*if(drawEllipses) {
		createEllipse(CALellipses[t]);
		transformEllipse(CALellipses[t],Transform);
		translateEllipse(CALellipses[t],Translate);
	}*/
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

	/*if(drawEllipses) {
		createEllipse(CALellipses[TIME_STEPS+1]);
		Matrix<3,3> Transform = zeros<3,3>();
		Transform(0,0) = xr; Transform(1,1) = yr; Transform(2,2) = zr;
		Matrix<3,1> Translate;
		Translate[0] = xc; Translate[1] = yc; Translate[2] = zc;
		transformEllipse(CALellipses[TIME_STEPS+1],Transform);
		translateEllipse(CALellipses[TIME_STEPS+1],Translate);
	}*/
}

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
	//REAL dist = gjk_distance(&VrelPoint, NULL, &ConvexHull, NULL, wptVrel, wptHull, NULL, 0); // GJK returns distance^2;

	//std::cout << "Distance: " << dist << std::endl;
	//std::cout << "Velocity Point: " << wptVrel[0] << " " << wptVrel[1] << " " << wptVrel[2] << std::endl;
	//std::cout << "Hull Point: " << wptHull[0] << " " << wptHull[1] << " " << wptHull[2] << std::endl;
	
	normal[0] = ((wptVrel[0]-wptHull[0])/(distance));
	normal[1] = ((wptVrel[1]-wptHull[1])/(distance));
	normal[2] = ((wptVrel[2]-wptHull[2])/(distance));
}

void pointInHull(const std::vector<Matrix<3,1>>& reachablePoints, bool& insideHull, const State& x0, const State&x1) {
	double distance;
	Matrix<3,1> normal;
	double epsilon = 0.0001;
	run_gjk(x0,x1,reachablePoints, distance, normal);
	if(distance < epsilon && distance > -1*epsilon)
		insideHull = true;
	else
		insideHull = false;
}

// Find the convex hull of the reachable points
void convexHull(const std::vector< Matrix<3,1>>& reachablePoints, bool& insideHull, const State& x0, const State& x1, double& distance, Matrix<3,1>& normal) {
	// Write the reachable pioints to a file to input into qconvex.exe
	std::ofstream outputFile;
	outputFile.open("pointList.txt");
	outputFile << "3" << std::endl << reachablePoints.size() << std::endl;
	for (int i = 0; i < reachablePoints.size(); i++)
		outputFile << reachablePoints[i][0] << " " << reachablePoints[i][1] << " " << reachablePoints[i][2] << std::endl;
	outputFile.close();

		// Use GJK, and if distance is 0, then the point is inside the hull.  If the point is inside the hull, the iterative distance method
		//     must be used to find the shortest vector out of the obstacle.
		// If the point is outside the hull, GJK can be used to find the shortest distance to the hull.
	system("qconvex n TO \"Planes.txt\" < \"pointList.txt\"");
	system("qconvex Fv TO \"facetVertices.txt\" < \"pointList.txt\"");
	Matrix<3,1> relVel;
	relVel[0] = x0[3]-x1[3]; relVel[1] = x0[4]-x1[4]; relVel[2] = x0[5]-x1[5];  // Relative velocity term

	std::vector<Matrix<4,1>> planeEq;	// Vector of plane equations	
	Matrix<4,1> tempPlane;				// Temp plane equation for pushback of planeEq
	std::string dump, first, second, third, fourth;	// Strings for reading values from text file
	int count = 0;						// dummy counter to skip first line of text read
	std::ifstream inputFile;
	inputFile.open("Planes.txt");
	int dim, numFacets;
	inputFile >> dim;
	inputFile >> numFacets;
	for (int i = 0; i < numFacets; ++i) {
		inputFile >> tempPlane[0];
		inputFile >> tempPlane[1];
		inputFile >> tempPlane[2];
		inputFile >> tempPlane[3];
		planeEq.push_back(tempPlane);
	}
	/*while(inputFile.good()) {
		if(count < 2) {
			getline(inputFile,dump);
			count++;
		}
		else {
			getline(inputFile,dump);
			std::istringstream line(dump);
			getline(line, first, ' ');
			getline(line, second,' ');
			getline(line, third, ' ');
			getline(line, fourth, ' ');
			tempPlane[0] = atof(first.c_str());
			tempPlane[1] = atof(second.c_str());
			tempPlane[2] = atof(third.c_str());
			tempPlane[3] = atof(fourth.c_str());
			planeEq.push_back(tempPlane);
		}
	}
	planeEq.pop_back();*/
	inputFile.close();

	std::vector<Matrix<3,1>> pointOnPlane;
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
	/*while(inputFile.good()) {
		if(count == 0) {
			getline(inputFile, dump);
			count++;
		}
		else {
			getline(inputFile,dump);
			std::istringstream line(dump);
			getline(line,first,' ');
			//index = atoi(first.c_str());
			tempPoint[0] = reachablePoints[atoi(first.c_str())][0]; tempPoint[1] = reachablePoints[atoi(first.c_str())][1]; tempPoint[2] = reachablePoints[atoi(first.c_str())][2];
			pointOnPlane.push_back(tempPoint);
		}
	}
	pointOnPlane.pop_back();*/
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

static const float RVO_EPSILON = 0.00001f;

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
        tLeft = std::max(tLeft, t);
      } else {
        /* Plane i bounds line on the right. */
        tRight = std::min(tRight, t);
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

  size_t linearProgram3(const std::vector<Plane>& planes, double radius, const Vector3& optVelocity, bool directionOpt, Vector3& result)
  {
    if (directionOpt) {
      /* 
       * Optimize direction. Note that the optimization velocity is of unit
       * length in this case.
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

  void computeHalfPlanes(const State& x0, const Matrix<3,1>& goalVel, Matrix<3,1>& newVel, const std::vector<double>& distance, const std::vector<Matrix<3,1>>& normals, const bool& insideHull){
	  
	  std::vector<Plane> orcaPlanes_; // Makes the half-plane vector of planes
	  Plane tempPlane;				// Used for filling the half-plane vector of planes
	  //for(int i = 0; i< distances.size(); i++) {
	  for(int i = 0; i < NUM_QUADS-1; i++) {
		  Vector3 tempNormal(normals[i][0], normals[i][1], normals[i][2]);
		  int mult;
		  if(insideHull)
			  //mult = 1;
			  mult = 1.01;
		  else
			  mult = -0.99;

		Vector3 tempVector(x0[3]+mult*distance[i]*normals[i][0], x0[4]+mult*distance[i]*normals[i][1], x0[5]+mult*distance[i]*normals[i][2]); 
		tempPlane.normal = tempNormal;
		tempPlane.point = tempVector;
		orcaPlanes_.push_back(tempPlane);
	  }

	  double maxSpeed_ = 10;
	  Vector3 prefVelocity_(goalVel[0],goalVel[1],goalVel[2]);
	  Vector3 newVelocity_;

	  size_t planeFail = linearProgram3(orcaPlanes_, maxSpeed_, prefVelocity_, false, newVelocity_);
    if (planeFail < orcaPlanes_.size()) {
      linearProgram4(orcaPlanes_, planeFail, maxSpeed_, newVelocity_);
    }

	newVel[0] = newVelocity_.x(); newVel[1] = newVelocity_.y(); newVel[2] = newVelocity_.z();
  }

int _tmain(int argc, _TCHAR* argv[]) {

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,U_DIM> B;
	Matrix<X_DIM,1> c;

	setup();
	std::vector< Quadrotor* > qlist(NUM_QUADS);

	Input uGoal;
	uGoal[0] = uGoal[1] = uGoal[2] = uGoal[3] = nominalInput;

	Qx = zeros<X_DIM,X_DIM>();

	Qv = zeros<V_DIM,V_DIM>();
	Qv(0,0) = Qv(1,1) = Qv(2,2) = 20;

	Qp = zeros<V_DIM,V_DIM>();
	Qp(0,0) = Qp(1,1) = Qp(2,2) = 20;

	R = 20*identity<U_DIM>();
	Rp = 20*identity<U_DIM>();

	M = 0.000000001 * identity<X_DIM>(); // motion noise variance
	N = 0.000000001 * identity<Z_DIM>(); // observation noise variance

	for( int i = 0; i < NUM_QUADS; i++) {
		qlist[i] = new Quadrotor();
		qlist[i]->setupQuadVisualization();
	}
		
	State xInit = zeros<X_DIM>();
	Rotation RotInit = identity<3>();
	Matrix<X_DIM,X_DIM> Pinit = 0.000000001 * identity<X_DIM>();
	State xGoal = zeros<X_DIM>();
	xGoal[12] = xGoal[13] = xGoal[14] = xGoal[15] = nominalInput;
	Rotation RotGoal = identity  <3>();
	//std::vector<Velocity> vGoalInit;
	Velocity VGoal0, VGoal1, VGoal2, VGoal3;
	Position PGoal0, PGoal1, PGoal2, PGoal3;
	xInit[12] = xInit[13] = xInit[14] = xInit[15] = nominalInput;

	// Quadrotor 0
	xInit[0] = -4.0; xInit[1] = 0.0; xInit[2] = 2.0;
	xInit[3] = 3.0; xInit[4] = 0.0; xInit[5] = 0;
	PGoal0[0]= 4.0; VGoal0[1] = 0.0; VGoal0[2] = 2.0;
	qlist[0]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, PGoal0, uGoal, RotGoal);
	// Quadrotor 1
	xInit[0] = 4.0; xInit[1] = 0.0; xInit[2] = 2.0;
	xInit[3] = -3.0; xInit[4] = 0.0; xInit[5] = 0;
	PGoal1[0] = 4.0; PGoal1[1] = 0.0; PGoal1[2] = 2.0;
	qlist[1]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, PGoal1, uGoal, RotGoal);
	/*// Quadrotor 2
	xInit[0] = 2.0; xInit[1] = 4.0; xInit[2] = 2.0;
	xInit[3] = -0.0; xInit[4] = -1.0; xInit[5] = 0.0;
	VGoal2[0] = -0.0; VGoal2[1] = -1.0; VGoal2[2] = 0.0;
	qlist[2]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, VGoal2, uGoal, RotGoal);
	// Quadrotor 3
	xInit[0] = -2.0; xInit[1] = -4.0; xInit[2] = 2;
	xInit[3] = 0.0; xInit[4] = 1.0; xInit[5] = 0;
	VGoal3[0] = 0.0; VGoal3[1] = 1.0; VGoal3[2] = 0;
	qlist[3]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, VGoal3, uGoal, RotGoal);*/
	/*double arm = sqrt(18.0);
	// Quadrotor 5
	xInit[0] = arm; xInit[1] = arm; xInit[2] = 2;
	xInit[3] = -arm/3.0; xInit[4] = -arm/3.0; xInit[5] = 0;
	tempVGoal[0] = -arm/3.0; tempVGoal[1] = -arm/3.0; tempVGoal[2] = 0;
	vGoalInit.push_back(tempVGoal);
	qlist[4]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, vGoalInit[4], uGoal, RotGoal);
	// Quadrotor 6
	xInit[0] = -arm; xInit[1] = arm; xInit[2] = 2;
	xInit[3] = arm/3.0; xInit[4] = -arm/3.0; xInit[5] = 0;
	tempVGoal[0] = arm/3.0; tempVGoal[1] = -arm/3.0; tempVGoal[2] = 0;
	vGoalInit.push_back(tempVGoal);
	qlist[5]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, vGoalInit[5], uGoal, RotGoal);
	// Quadrotor 7
	xInit[0] = -a; xInit[1] = -a; xInit[2] = 2;
	xInit[3] = a/3.0; xInit[4] = a/3.0; xInit[5] = 0;
	tempVGoal[0] = a/3.0; tempVGoal[1] = a/3.0; tempVGoal[2] = 0;
	vGoalInit.push_back(tempVGoal);
	qlist[6]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, vGoalInit[6], uGoal, RotGoal);
	//Quadrotor 8
	xInit[0] = a; xInit[1] = -a; xInit[2] = 2;
	xInit[3] = -a/3.0; xInit[4] = a/3.0; xInit[5] = 0;
	tempVGoal[0] = -a/3.0; tempVGoal[1] = a/3.0; tempVGoal[2] = 0;
	vGoalInit.push_back(tempVGoal);
	qlist[7]->setupQuadrotors(xInit, RotInit, Pinit, xGoal, vGoalInit[7], uGoal, RotGoal);*/

	Matrix<3,X_DIM> C = zeros<3,X_DIM>();			// Maps state to position
	C(0,0) = C(1,1) = C(2,2) = 1;
	Matrix<V_DIM,X_DIM> V = zeros<V_DIM,X_DIM>();   // Maps state to velocity
	V(0,3) = V(1,4) = V(2,5) = 1;

	Matrix<X_DIM,X_DIM> Ft = identity<X_DIM>();		// Ft = A^t
	Matrix<X_DIM,V_DIM> Gt = zeros<X_DIM,V_DIM>();  // Gt = (A^(t-1) + A^(t-2) + A^(t-3) ... + A + 1)B

	CAL_DestroyGroup(rotordummy);
	CAL_DestroyGroup(ellipsedummy);

	for (size_t i = 0; i < NUM_QUADS; ++i) {
		qlist[i]->findMatrices(A,B,c);
		qlist[i]->vGoal = qlist[i]->findVGoal();
	}

	std::vector< Matrix<3,1> > points;		// Vector to save the points of a sampled sphere
	createSpheres(points);					// Function that samples the sphere
	std::vector< Matrix<3,1>> hullVertices;	// Vector of the vertices of the convex hull
	bool insideHull = false;				// Bool to test if the relative velocity point is inside or outside the convex hull
	int numVertices;						// Number of vertices of the convex hull
	Matrix<4,1> halfPlaneVector;         // Elements 0-2 = normal vector out of the convex hull, element 3 = magnitude of distance vector
	double distance;					// Distance for halfPlaneVector;
	std::vector<double> distances;
	Matrix<3,1> normalVector;			// Normal vector for halfPlaneVector
	std::vector<Matrix<3,1>> normalVectors;
	double maxSpeed = 30.0;
	//int counter = 0;

	int obj = 0;
	for(size_t t = 0; t < TIME_STEPS; ++t) {
		std::cout << t << std::endl;
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 0 to 1
		distances.clear();
		normalVectors.clear();
		qlist[0]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[0]->L, qlist[0]->E, Ft, Gt);
			createObstacle(qlist[0]->reachablePoints,points,qlist[0]->Transform,qlist[0]->Translate,t,Ft,Gt,C,qlist[0]->x,qlist[1]->x);
		}
		pointInHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[1]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[1]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[0]->x, qlist[1]->x, qlist[0]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 0 to 2
		qlist[0]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[0]->L, qlist[0]->E, Ft, Gt);
			createObstacle(qlist[0]->reachablePoints,points,qlist[0]->Transform,qlist[0]->Translate,t,Ft,Gt,C,qlist[0]->x,qlist[2]->x);
		}
		pointInHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[2]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[2]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[0]->x, qlist[2]->x, qlist[0]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);*/
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 0 to 3
		qlist[0]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[0]->L, qlist[0]->E, Ft, Gt);
			createObstacle(qlist[0]->reachablePoints,points,qlist[0]->Transform,qlist[0]->Translate,t,Ft,Gt,C,qlist[0]->x,qlist[3]->x);
		}
		pointInHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[3]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[0]->reachablePoints, insideHull, qlist[0]->x, qlist[3]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[0]->x, qlist[3]->x, qlist[0]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);*/
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 0 Half Planes
		computeHalfPlanes(qlist[0]->x, qlist[0]->vGoal, qlist[0]->newV, distances, normalVectors, insideHull);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 1 to 0
		distances.clear();
		normalVectors.clear();
		qlist[1]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[1]->L, qlist[1]->E, Ft, Gt);
			createObstacle(qlist[1]->reachablePoints,points,qlist[1]->Transform,qlist[1]->Translate,t,Ft,Gt,C,qlist[1]->x,qlist[0]->x);
		}
		pointInHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[0]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[0]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[1]->x, qlist[0]->x, qlist[1]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 1 to 2
		qlist[1]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[1]->L, qlist[1]->E, Ft, Gt);
			createObstacle(qlist[1]->reachablePoints,points,qlist[1]->Transform,qlist[1]->Translate,t,Ft,Gt,C,qlist[1]->x,qlist[2]->x);
		}
		pointInHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[2]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[2]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[1]->x, qlist[2]->x, qlist[1]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);*/
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 1 to 3
		qlist[1]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[1]->L, qlist[1]->E, Ft, Gt);
			createObstacle(qlist[1]->reachablePoints,points,qlist[1]->Transform,qlist[1]->Translate,t,Ft,Gt,C,qlist[1]->x,qlist[3]->x);
		}
		pointInHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[3]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[1]->reachablePoints, insideHull, qlist[1]->x, qlist[3]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[1]->x, qlist[3]->x, qlist[1]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);*/
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 1 Half Planes
		computeHalfPlanes(qlist[1]->x, qlist[1]->vGoal, qlist[1]->newV, distances, normalVectors, insideHull);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 2 to 0
		distances.clear();
		normalVectors.clear();
		qlist[2]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[2]->L, qlist[2]->E, Ft, Gt);
			createObstacle(qlist[2]->reachablePoints,points,qlist[2]->Transform,qlist[2]->Translate,t,Ft,Gt,C,qlist[2]->x,qlist[0]->x);
		}
		pointInHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[0]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[0]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[2]->x, qlist[0]->x, qlist[2]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 2 to 1
		qlist[2]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[2]->L, qlist[2]->E, Ft, Gt);
			createObstacle(qlist[2]->reachablePoints,points,qlist[2]->Transform,qlist[2]->Translate,t,Ft,Gt,C,qlist[2]->x,qlist[1]->x);
		}
		pointInHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[1]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[1]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[2]->x, qlist[1]->x, qlist[2]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 2 to 3
		qlist[2]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[2]->L, qlist[2]->E, Ft, Gt);
			createObstacle(qlist[2]->reachablePoints,points,qlist[2]->Transform,qlist[2]->Translate,t,Ft,Gt,C,qlist[2]->x,qlist[3]->x);
		}
		pointInHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[3]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[2]->reachablePoints, insideHull, qlist[2]->x, qlist[3]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[2]->x, qlist[3]->x, qlist[2]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 2 Half Planes
		computeHalfPlanes(qlist[2]->x, qlist[2]->vGoal, qlist[2]->newV, distances, normalVectors, insideHull);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 3 to 0
		distances.clear();
		normalVectors.clear();
		qlist[3]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[3]->L, qlist[3]->E, Ft, Gt);
			createObstacle(qlist[3]->reachablePoints,points,qlist[3]->Transform,qlist[3]->Translate,t,Ft,Gt,C,qlist[3]->x,qlist[0]->x);
		}
		pointInHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[0]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[0]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[3]->x, qlist[0]->x, qlist[3]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 3 to 1
		qlist[3]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[3]->L, qlist[3]->E, Ft, Gt);
			createObstacle(qlist[3]->reachablePoints,points,qlist[3]->Transform,qlist[3]->Translate,t,Ft,Gt,C,qlist[3]->x,qlist[1]->x);
		}
		pointInHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[1]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[1]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[3]->x, qlist[1]->x, qlist[3]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 3 to 2
		qlist[3]->reachablePoints.clear();
		Ft = identity<X_DIM>();
		Gt = zeros<X_DIM,V_DIM>();
		for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
			findFG(A,B,qlist[3]->L, qlist[3]->E, Ft, Gt);
			createObstacle(qlist[3]->reachablePoints,points,qlist[3]->Transform,qlist[3]->Translate,t,Ft,Gt,C,qlist[3]->x,qlist[2]->x);
		}
		pointInHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[2]->x);
		if(insideHull) { // if the point is inside the convex hull
			convexHull(qlist[3]->reachablePoints, insideHull, qlist[3]->x, qlist[2]->x, distance, normalVector);
		}
		else // If the point is outside the convex hull, run GJK
		{
			run_gjk(qlist[3]->x, qlist[2]->x, qlist[3]->reachablePoints, distance, normalVector);
		}	
		distances.push_back(0.5*distance);
		normalVectors.push_back(normalVector);
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 3 Half Planes
		computeHalfPlanes(qlist[3]->x, qlist[3]->vGoal, qlist[3]->newV, distances, normalVectors, insideHull);*/
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 0 Update
		qlist[0]->vGoal = qlist[0]->newV;
		Input u = qlist[0]->findU();
		qlist[0]->propagateU(u);
		kalmanFilter1(qlist[0]->x, qlist[0]->Rot, u, M, qlist[0]->P);
		Observation z = sampleGaussian(h(qlist[0]->xTrue, qlist[0]->RotTrue), N);
		kalmanFilter2(qlist[0]->x, qlist[0]->Rot, z, N, qlist[0]->P);
		qlist[0]->visualize(t*dt);
		qlist[0]->vGoal = qlist[0]->findVGoal();
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 1 Update
		qlist[1]->vGoal = qlist[1]->newV;
		u = qlist[1]->findU();
		qlist[1]->propagateU(u);
		kalmanFilter1(qlist[1]->x, qlist[1]->Rot, u, M, qlist[1]->P);
		z = sampleGaussian(h(qlist[1]->xTrue, qlist[1]->RotTrue), N);
		kalmanFilter2(qlist[1]->x, qlist[1]->Rot, z, N, qlist[1]->P);
		qlist[1]->visualize(t*dt);
		qlist[1]->vGoal = qlist[1]->findVGoal();
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		/*// Rotor 2 Update
		qlist[2]->vGoal = qlist[2]->newV;
		u = qlist[2]->findU();
		qlist[2]->propagateU(u);
		kalmanFilter1(qlist[2]->x, qlist[2]->Rot, u, M, qlist[2]->P);
		z = sampleGaussian(h(qlist[2]->xTrue, qlist[2]->RotTrue), N);
		kalmanFilter2(qlist[2]->x, qlist[2]->Rot, z, N, qlist[2]->P);
		qlist[2]->visualize(t*dt);
		qlist[2]->vGoal = VGoal2;
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------
		// Rotor 3 Update
		qlist[3]->vGoal = qlist[3]->newV;
		u = qlist[3]->findU();
		qlist[3]->propagateU(u);
		kalmanFilter1(qlist[3]->x, qlist[3]->Rot, u, M, qlist[3]->P);
		z = sampleGaussian(h(qlist[3]->xTrue, qlist[3]->RotTrue), N);
		kalmanFilter2(qlist[3]->x, qlist[3]->Rot, z, N, qlist[3]->P);
		qlist[3]->visualize(t*dt);
		qlist[3]->vGoal = VGoal3;*/
	}

	/*std::ofstream outputFile;
	outputFile.open("mainTest.txt");
	for(size_t t = 0; t < TIME_STEPS; ++t) {	// Step through dt (dt = 1/30s )
		//for(int i = 0; i < 1; ++i)
			//std::cout << "Rotor [" << i << "] Initial VGoal: " << vGoalInit[i][0] << " " << vGoalInit[i][1] << " " << vGoalInit[i][2] << std::endl;
		for( size_t i = 0; i< NUM_QUADS; ++i) {	// Step through all the quadrotors to update their motion and run collision avoidance
			start = clock();
			normalVector = zeros<3,1>();
			distance = 0;
			distances.clear();
			normalVectors.clear();
			//fstart = clock();
			//qlist[i]->ellipsoids.clear();
			//qlist[i]->reachablePoints.clear();
			for(size_t j = 0; j < NUM_QUADS; ++j) { // Loop through all the quadrotors and make an obstacle between quadrotor i and quadrotor j to update the velocity
				if(i==j);
				else {
				// Clears ellipsoid and reachable point vectors to prepare for next obstacle creation
					qlist[i]->ellipsoids.clear();
					qlist[i]->reachablePoints.clear();
				// Loop through all obstacle time steps and make the obstacle ellipses
					Ft = identity<X_DIM>();
					Gt = zeros<X_DIM,V_DIM>();
					for(size_t k = 0; k < OBSTACLE_STEPS; ++k) {
						findFG(A,B,qlist[i]->L, qlist[i]->E, Ft, Gt);
						createObstacle(qlist[i]->reachablePoints,points,qlist[i]->Transform,qlist[i]->Translate,t,Ft,Gt,C,qlist[i]->x,qlist[j]->x,1);
					}
				// Once those are made, find the reachable points in the obstacle ellipses
					//findReachableObstacle(qlist[i]->ellipsoids,qlist[i]->reachablePoints,qlist[i]->x,qlist[j]->x, maxSpeed);
					//qlist[i]->reachablePoints = qlist[i]->ellipsoids;
				// If the reachable points are enough to make a hull, run convex hull
					//if(qlist[i]->reachablePoints.size() > 4) {
					pointInHull(qlist[i]->reachablePoints, insideHull, qlist[i]->x, qlist[j]->x);
					if(insideHull) { // if the point is inside the convex hull
						convexHull(qlist[i]->reachablePoints, insideHull, qlist[i]->x, qlist[j]->x, distance, normalVector);
					}
					else // If the point is outside the convex hull, run GJK
					{
						run_gjk(qlist[i]->x, qlist[j]->x, qlist[i]->reachablePoints, distance, normalVector);
					}	
					
					//std::cout << "t: " << t << ", distance = " << distance << ", normal: " << normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] <<std::endl;
					/*outputFile << "Quadrotor: " << i << std::endl;
					outputFile << "t: " << t << std::endl << "distance: "<< distance << std::endl << "normal: " << normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << std::endl;
					outputFile << "insideHull: " << insideHull << std::endl;
					outputFile << std::endl << std::endl;
					
					distances.push_back(0.5*distance);
					normalVectors.push_back(normalVector);
					if(i ==20) {
						std::cout << "Rotor[" << i <<"] Distance: " << distances[0] << std::endl;
						std::cout << "Rotor[" << i <<"] Normal V: " << normalVectors[0][0] << " " << normalVectors[0][1] << " " << normalVectors[0][2] << std::endl;
						std::cout << "Rotor[" << i <<"] insideHu: " << insideHull << std::endl << std::endl;
					}
					//    create the half-plane and find a new velocity that is within the safe zone of the half planes
					//    Save the new velocity as a separate variable, and keep the old one
					//    move on to the next quadrotor using the old velocity to make it's obstacle, once all quadrotors have been
					//        collision checked, it is safe to update all of their velocities	
				}
			}
			computeHalfPlanes(qlist[i]->x, qlist[i]->vGoal, qlist[i]->newV, distances, normalVectors, insideHull);
			//if(i==0)
				//std::cout << "Quadrotor " << i << "new goal v: " << qlist[i]->newV[0] << " " << qlist[i]->newV[1] << " " << qlist[i]->newV[2] <<std::endl;
		}

		for(int i = 0; i < NUM_QUADS; i++) {
				//qlist[i]->setupQuadrotors(qlist[i]->x, qlist[i]->Rot, qlist[i]->P, qlist[i]->xGoal, qlist[i]->newV, qlist[i]->uGoal, qlist[i]->RotGoal); 
				qlist[i]->vGoal = qlist[i]->newV;
				Input u = qlist[i]->findU();
			// Send control input to quadrotor
				qlist[i]->propagateU(u);
			// Process control input in time update of Kalman filter
				kalmanFilter1(qlist[i]->x, qlist[i]->Rot, u, M, qlist[i]->P);
			// Receive observation
				Observation z = sampleGaussian(h(qlist[i]->xTrue, qlist[i]->RotTrue), N);
			// Process observation in measurement of Kalman filter
				kalmanFilter2(qlist[i]->x, qlist[i]->Rot, z, N, qlist[i]->P);
			// Find the F(t) and G(t) used to create the obstacles
				//findFG(A,B,qlist[i]->L, qlist[i]->E, Ft, Gt);
				qlist[i]->visualize(t*dt);
				qlist[i]->vGoal = vGoalInit[i];
			/*else if(i==3)
				qlist[i]->vGoal = tempVGoal3;
			else
				qlist[i]->vGoal = tempVGoal4;
		}
		for(int i = 0; i<NUM_QUADS; i++) qlist[i]->vGoal = vGoalInit[i];
		stop = clock();
		std::cout << t <<std::endl;
		//std::cout << "Duration of time step for one quadrotor " << t << ": " << (stop-start) - (fstop-fstart) << std::endl << std::endl;
		// Repeat process to update the collision-free velocity of all quadrotors
		
	} // Iterate that process through the full simulation
	*/
	std::cout << "Done" << std::endl;
	int k;
	std::cin >> k;

	CAL_End();
	return 0;
}
