#ifndef TIMEINTEGRATOR_H
#define TIMEINTEGRATOR_H

#include "solver.h"

namespace PDtools {
//------------------------------------------------------------------------------
class TimeIntegrator : public Solver {
protected:
  int m_indexRho;

public:
  TimeIntegrator() { ; }
  virtual ~TimeIntegrator() { ; }

  virtual void solve();
  virtual void stepForward(int i);
  enum IntegratorErrorMessages { TimeStepNotSet, ParticlesMissingRhoOrMass };

protected:
  virtual void checkInitialization();
  virtual void initialize();
  virtual void integrateStepOne() = 0;
  virtual void integrateStepTwo() = 0;
};
//------------------------------------------------------------------------------
}
#endif // TIMEINTEGRATOR_H
