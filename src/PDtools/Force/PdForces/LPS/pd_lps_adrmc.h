#ifndef PD_LPS_ADRMC_H
#define PD_LPS_ADRMC_H

#include "PDtools/Force/PdForces/LPS/lps_mc.h"

namespace PDtools {
//------------------------------------------------------------------------------
class PD_LPS_adrmc : public LPS_mc {
public:
  PD_LPS_adrmc(PD_Particles &particles, double phi, double C, double T,
               bool planeStress, bool analyticalM = false);

  virtual void calculateForces(const int id, const int i);
  virtual void evaluateStatic(int id, int i);
  virtual void evaluateStepTwo(int id, int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_ADRMC_H
