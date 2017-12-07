#ifndef PD_LPS_ADRMC_POROSITY_P_H
#define PD_LPS_ADRMC_POROSITY_P_H

#include "PDtools/Force/PdForces/LPS_porosity/lps_p_mc.h"

namespace PDtools {
//------------------------------------------------------------------------------
class PD_LPS_porosity_adrmc : public LPS_porosity_mc {
public:
  PD_LPS_porosity_adrmc(PD_Particles &particles, double m, double b, double phi,
                        double C, double T, bool planeStress,
                        bool analyticalM = false);

  virtual void calculateForces(const int id, const int i);

  virtual void evaluateStatic(int id, int i);

  virtual void evaluateStepTwo(int id, int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_ADRMC_POROSITY_P_H
