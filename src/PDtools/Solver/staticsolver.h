#ifndef STATICSOLVER_H
#define STATICSOLVER_H

#include "solver.h"

namespace PDtools {
//------------------------------------------------------------------------------
class StaticSolver : public Solver {
public:
  StaticSolver(int maxIterations, int threshold);

  virtual void solve();
  virtual void stepForward(int i);
  virtual void iterate();

protected:
  int m_maxIterations;
  double m_threshold;
  int m_degFreedom;
  arma::sp_mat C;
  vec u_k;
  vec r_k;
  vec r_k1;
  vec p_k;
  vec Ap;
  vec b_k;

  virtual void initialize();
  virtual void save(int i);
  void createStiffnessMatrix();
  void updateStiffnessMatrix();
  void computeStress();
};
//------------------------------------------------------------------------------
}
#endif // STATICSOLVER_H
