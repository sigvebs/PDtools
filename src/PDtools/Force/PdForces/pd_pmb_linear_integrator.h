#ifndef PD_PMB_LINEAR_INTEGRATOR_H
#define PD_PMB_LINEAR_INTEGRATOR_H

#include "PDtools/Force/PdForces/pd_bondforce.h"

//------------------------------------------------------------------------------
namespace PDtools {

class PD_PMB_LINEAR_INTEGRATOR : public PD_bondForce {
public:
  PD_PMB_LINEAR_INTEGRATOR(PD_Particles &particles, double lc, double delta,
                           double alpha);

  virtual void calculateForces(const int id, const int i);

protected:
  double m_lc;
  double m_delta;
  double m_alpha;
};
//------------------------------------------------------------------------------
}

#endif // PD_PMB_LINEAR_INTEGRATOR_H
