#ifndef VELOCITYVERLETINTEGRATOR_H
#define VELOCITYVERLETINTEGRATOR_H

#include "PDtools/Solver/timeintegrator.h"

namespace PDtools {
//------------------------------------------------------------------------------
class VelocityVerletIntegrator : public TimeIntegrator {
public:
  VelocityVerletIntegrator();
  ~VelocityVerletIntegrator();

  virtual void integrateStepOne();
  virtual void integrateStepTwo();
};
//------------------------------------------------------------------------------
}
#endif // VELOCITYVERLETINTEGRATOR_H
