#include "simulator.h"

#define NUM_QUADS 1

/*#define X_DIM 16         // State dimension
#define V_DIM 3
#define U_DIM 4          // Control input dimension
#define Z_DIM 6          // Observation dimension

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
int rotordummy;
int ellipsedummy;
int cal_scene;
int cal_floor;
int cal_ellipse;

Matrix<X_DIM,X_DIM> Qx;
Matrix<V_DIM,V_DIM> Qv;
Matrix<U_DIM,U_DIM> R;
Matrix<X_DIM,X_DIM> P;
Matrix<X_DIM,X_DIM> M;
Matrix<Z_DIM,Z_DIM> N;

// Quadrotor Class
class Quadrotor {
		//State xInit;
		State xGoal;
		Velocity vGoal;
		Input uGoal;
		Rotation RotGoal;
		//Matrix<X_DIM,X_DIM> Qx;
		//Matrix<V_DIM,V_DIM> Qv;
		//Matrix<U_DIM,U_DIM> R;
		//Matrix<X_DIM,X_DIM> P;
		//Matrix<X_DIM,X_DIM> M;
		//Matrix<Z_DIM,Z_DIM> N;
		int cal_quadrotor;
		int cal_ellipse;

		State x;
	    Rotation Rot;
	    State xTrue;
	    Rotation RotTrue;
	    Matrix<U_DIM,X_DIM> L;
	    Matrix<U_DIM,V_DIM> E;
	    Matrix<U_DIM,1> l; // ell*/

	public:
		void setupQuadrotors(const State& xfinal, const Velocity& vfinal, const Input& uHover, const Rotation& Rgoal) {
			xGoal = xfinal;
			vGoal = vfinal;
			uGoal = uHover;
			RotGoal = Rgoal;
			//P = Pinit;
		}

		void setupQuadVisualization() {
			CAL_CloneGroup(&cal_quadrotor, rotordummy, 0, false, "Quadrotor");
			CAL_CloneGroup(&cal_ellipse, ellipsedummy, 0, false, "ellipse");
		}

		void findMatrices(const Rotation& RotInit) {
			Rot = RotInit;
			xTrue = sampleGaussian(x,P);
			RotTrue = RotInit;
			RotTrue = RotTrue * exp(skewSymmetric((xTrue).subMatrix<3,1>(6,0)));
			xTrue[6] = 0; xTrue[7] = 0; xTrue[8] = 0;
			controlMatrices(uGoal, xGoal, L, E, l);
		}

void visualize(double t) {
    float p[3] = {(float) xTrue[0], (float) xTrue[1], (float) xTrue[2]};
    Matrix<4> q = quatFromRot(RotTrue);
    float o[4] = {(float) q[0], (float) q[1], (float) q[2], (float) q[3]};
    CAL_AddGroupKeyState(cal_quadrotor, (float) t, p, o);
    //CAL_SetGroupPosition(cal_quadrotor, p[0], p[1], p[2]);
    //CAL_SetGroupQuaternion(cal_quadrotor, q[0], q[1], q[2], q[3]);

    std::cout << ~(xTrue.subMatrix<3>(0,3));

    /*Matrix<3,3> V, E;
    jacobi(P.subMatrix<3,3>(0,0), V, E);
    q = quatFromRot(V);
    float pos[3] = {(float) xEst[0], (float) xEst[1], (float) xEst[2]};
    float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
    float scale[3] = {(float) (2*sqrt(E(0,0))), (float) (2*sqrt(E(1,1))), (float) (2*sqrt(E(2,2)))};
    CAL_AddGroupKeyState(cal_ellipse, (float) t, pos, quat, scale);*/

    Matrix<3,3> V, E;
    jacobi(P.subMatrix<3,3>(0,0), V, E);
    q = quatFromRot(V);
    float pos[3] = {(float) x[0], (float) x[1], (float) x[2]};
    float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
    float scale[3] = {(float) (2*sqrt(E(0,0))), (float) (2*sqrt(E(1,1))), (float) (2*sqrt(E(2,2)))};
    CAL_AddGroupKeyState(cal_ellipse, (float) t, pos, quat, scale);
  }
	
		Input findU() {
			return riccatiControllerSteady(x, Rot, vGoal, RotGoal, uGoal, L, E, l);
		}
		
		void propagateU( const Input& u ) {
			propagate(xTrue, RotTrue, u, M);
		}

		State getState() {
			return x;
		}

		State getTrueState() {
			return xTrue;
		}

		Rotation getRot() {
			return Rot;
		}

		Rotation getTrueRot() {
			return RotTrue;
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

  double volumeWidth   = 6.8834; // m
  double volumeLength  = 4.572;  // m
  double volumeHeight  = 3.3528;  // m

  // visualization
  CAL_Initialisation();
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
  CAL_CreateSphere(ellipsedummy, 1, 0, 0, 0);
  
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

inline void controlMatrices(const Input& uGoal, const Matrix<X_DIM,1>& xstar, Matrix<U_DIM,X_DIM>& L, Matrix<U_DIM,V_DIM>& E, Matrix<U_DIM,1>& l) {
	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,U_DIM> B;
	Matrix<X_DIM,1> c;
  
	Matrix<V_DIM,X_DIM> V = zeros<V_DIM,X_DIM>();
	V(0,3) = V(1,4) = V(2,5) = 1;

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
	Matrix<X_DIM, X_DIM> S = ~V*Qv*V + Qx;
	for (int i = 0; i < 200; ++i) {
		//S = Q + ~A*S*A - ~A*S*B*!(R + ~B*S*B)*~B*S*A;
		S = ~V*Qv*V + Qx + ~A*S*A - ~A*S*B*!(R + ~B*S*B)*(~B*S*A);
	}

	//Matrix<X_DIM,X_DIM> T = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * Q;
    //Matrix<X_DIM,V_DIM> T = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>())*(~V*Qv);
	Matrix<X_DIM,V_DIM> T = -S*~V;
	//Matrix<X_DIM,1>     v = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>()) * (~A*S*B*!(R + ~B*S*B)*~B*S - ~A*S)*c;
    Matrix<X_DIM,1>     a = pseudoInverse(~A - ~A*S*B*!(R + ~B*S*B)*~B - identity<X_DIM>())*(Qx*xstar - ~A*S*c + ~A*S*B*!(R+~B*S*B)*~B*S*c);

	// feedback law: u = L*xTilde + E*xTildeGoal + l
	L = -!(R + ~B*S*B)*~B*S*A;
	E = -!(R + ~B*S*B)*~B*T;
	l = -!(R + ~B*S*B)*(~B*S*c + ~B*a); // should be zero for steady-state linearization
}

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
  
  /*State xTildeGoal; //= xGoal - xHat;
  xTildeGoal.insert(0,0, ~RLocal * xGoal.subMatrix<3,1>(0,0));
  xTildeGoal.insert(3,0, ~RLocal * xGoal.subMatrix<3,1>(3,0));
  xTildeGoal.insert(6,0, errFromRot(~RLocal*RGoal));
  xTildeGoal.insert(9,0, xGoal.subMatrix<3,1>(9,0));
  xTildeGoal.insert(12,0, xGoal.subMatrix<4,1>(12,0) - uGoal);*/

  Velocity vTildeGoal;
  vTildeGoal.insert(0,0, ~RLocal * vGoal);

  State xTilde; //= x - xHat;
  xTilde.insert(0,0, ~RLocal * x.subMatrix<3,1>(0,0));
  xTilde.insert(3,0, ~RLocal * x.subMatrix<3,1>(3,0));
  xTilde.insert(6,0, errFromRot(~RLocal*R0)); //xTilde.insert(6,0, -(~R0*axis));
  xTilde.insert(9,0, x.subMatrix<3,1>(9,0));
  xTilde.insert(12,0, x.subMatrix<4,1>(12,0) - uGoal);
  
  return uGoal + L*xTilde + E*vTildeGoal + l;
}

/*inline Input riccatiControllerCurrent(const State& x, const Rotation& R0, const State& xGoal, const Rotation& RGoal, const Input& uGoal, const Matrix<X_DIM, X_DIM>& Q, const Matrix<U_DIM, U_DIM>& R, Matrix<X_DIM, X_DIM>& S) {
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
}*/

int _tmain(int argc, _TCHAR* argv[])
{
	setup();
	std::vector< Quadrotor* > qlist(NUM_QUADS);

	Input uGoal;
	uGoal[0] = uGoal[1] = uGoal[2] = uGoal[3] = nominalInput;

	Qx = zeros<X_DIM,X_DIM>();
  // Position & Velocity = 0;
	//Qx(6,6) = Qx(7,7) = 1; Qx(8,8) = 20;
	//Qx(9,9) = Qx(10,10) = Qx(11,11) = 1;

	Qv = zeros<V_DIM,V_DIM>();
	Qv(0,0) = Qv(1,1) = Qv(2,2) = 5;

	R = 20*identity<U_DIM>();

	M = 0.0001 * identity<X_DIM>(); // motion noise variance
	N = 0.0001 * identity<Z_DIM>(); // observation noise variance
	P = 0.001 * identity<X_DIM>();

	uGoal[0] = uGoal[1] = uGoal[2] = uGoal[3] = nominalInput;
	
	for( int i = 0; i < NUM_QUADS; i++) {
		qlist[i] = new Quadrotor();
		qlist[i]->setupQuadVisualization();

		State xInit = zeros<X_DIM>(); 
		//xInit[0] = random()*6 - 3; xInit[1] = random()*4 - 2;
		xInit[0] = 0; xInit[1] = 0; xInit[3] = 0;
		xInit[12] = xInit[13] = xInit[14] = xInit[15] = nominalInput;

		State xGoal = zeros<X_DIM>();
		xGoal[12] = xGoal[13] = xGoal[14] = xGoal[15] = nominalInput;

		Velocity vGoal = zeros<V_DIM,1>();
		vGoal[0] = 0; vGoal[1]= 0; vGoal[2] = 0;

		std::cout << i << std::endl << ~xInit << std::endl << ~xGoal << std::endl << std::endl; 

		qlist[i]->setupQuadrotors(xGoal, vGoal, uGoal, identity<3>());
	}

	// Destroys the dummy initializations of quadrotor and ellipse
	//     so that they do not appear at the origin
	CAL_DestroyGroup(rotordummy);
	CAL_DestroyGroup(ellipsedummy);

	clock_t beginTime = clock();

	for (size_t i = 0; i < NUM_QUADS; ++i) {
		qlist[i]->findMatrices(identity<3>());

		/*qlist[i]->x = qlist[i]->xInit;
		qlist[i]->Rot = qlist[i]->RotInit;
		qlist[i]->xTrue = sampleGaussian(qlist[i]->xInit, qlist[i]->P);
		qlist[i]->RotTrue = qlist[i]->RotInit;
		qlist[i]->RotTrue = qlist[i]->RotTrue*exp(skewSymmetric(qlist[i]->xTrue.subMatrix<3,1>(6,0)));
		qlist[i]->xTrue[6] = 0; qlist[i]->xTrue[7] = 0; qlist[i]->xTrue[8] = 0;

		controlMatrices(qlist[i]->uGoal, V, qlist[i]->Qv, qlist[i]->Qx, qlist[i]->xGoal, qlist[i]->R, qlist[i]->L, qlist[i]->E, qlist[i]->l);*/
	}

	  int k;
	  std::cin >> k;

    for (size_t t = 0; t < 320; ++t) {
		for (size_t i = 0; i < NUM_QUADS; ++i) {
			Input u = qlist[i]->findU();
   
			// Send control input to quadrotor
			qlist[i]->propagateU(u);
			//propagate(qlist[i]->xTrue, qlist[i]->RotTrue, u, M);
			// Process control input in time update of Kalman filter
			kalmanFilter1(qlist[i]->getState(), qlist[i]->getRot(), u, M, P);

			// Receive observation
			Observation z = sampleGaussian(h(qlist[i]->getTrueState(), qlist[i]->getTrueRot()), N);
			// Process observation in measurement of Kalman filter
			kalmanFilter2(qlist[i]->getState(), qlist[i]->getRot(), z, N, P);

			qlist[i]->visualize(t*dt);
		}
	}

  std::cout << (clock() - beginTime) / (double) (320*CLOCKS_PER_SEC) << std::endl;


  std::cin >> k;
  
  CAL_End();
 	return 0;
}