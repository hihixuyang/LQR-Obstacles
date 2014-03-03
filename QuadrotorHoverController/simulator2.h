#include "stdafx.h"
#include <vector>

#define X_DIM 16         // State dimension
#define V_DIM 3
#define U_DIM 4          // Control input dimension
#define Z_DIM 6          // Observation dimension

typedef Matrix<X_DIM> State;
typedef Matrix<U_DIM> Input;
typedef Matrix<Z_DIM> Observation;
typedef Matrix<3,3>   Rotation;
typedef Matrix<3,1>   Velocity;
typedef Matrix<3,1>   Position;

inline double random();
inline double normal();

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

inline State f(const State&, const Rotation&, const Input&);
inline Observation h(const State&, const Rotation&);
inline Matrix<X_DIM,X_DIM> Jacobian_fx(const State&, const Rotation&, const Input&);
inline Matrix<X_DIM, U_DIM> Jacobian_fu(const State&, const Rotation&, const Input&);
inline Matrix<Z_DIM, X_DIM> Jacobian_hx(const State&, const Rotation&);
inline void linearizeDiscretize(const State&, const Rotation&, const Input&, Matrix<X_DIM,X_DIM>&, Matrix<X_DIM,U_DIM>&, Matrix<X_DIM,1>&);
inline void propagate(State&, Rotation&, const Input&, const Matrix<X_DIM,X_DIM>&);
inline void kalmanFilter1(State&, Rotation&, const Input&, const Matrix<X_DIM,X_DIM>&, Matrix<X_DIM,X_DIM>&);
inline void kalmanFilter2(State&, Rotation&, const Observation&, const Matrix<Z_DIM, Z_DIM>&, Matrix<X_DIM, X_DIM>&);
inline void controlMatrices(const Input& uGoal, const Matrix<X_DIM,1>& xstar, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,1>& c, Matrix<U_DIM,X_DIM>& L, Matrix<U_DIM,V_DIM>& E, Matrix<U_DIM,1>& l, Matrix<V_DIM,X_DIM>& Lh, Matrix<V_DIM,V_DIM>& Eh);
inline Input riccatiControllerSteady(const State&, const Rotation&, const Velocity&, const Rotation&, const Input&, const Matrix<U_DIM,X_DIM>&, const Matrix<U_DIM,V_DIM>&, const Matrix<U_DIM,1>&);
inline Velocity riccatiControllerSteadyPosition(const State&, const Rotation&, const Position&, const Rotation&, const Input&, const Matrix<V_DIM,X_DIM>&, const Matrix<V_DIM,V_DIM>&);

