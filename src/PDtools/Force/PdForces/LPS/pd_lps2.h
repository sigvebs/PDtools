#ifndef PD_LPS2_H
#define PD_LPS2_H

#include "PDtools/Force/PdForces/LPS/pd_lps.h"
namespace PDtools {
//------------------------------------------------------------------------------
class PD_LPS2 : public PD_LPS {
public:
  PD_LPS2(PD_Particles &particles, bool planeStress = false,
          bool analyticalM = false);

  virtual void calculateForces(const int id, const int i);
  virtual void evaluateStepOne();
};
//------------------------------------------------------------------------------
}

#endif // PD_LPS2_H
