#ifndef PD_LPS_ADR_STRAIN_H
#define PD_LPS_ADR_STRAIN_H

#include "PDtools/Force/PdForces/LPS/pd_lps_crit_strain.h"

namespace PDtools {
//------------------------------------------------------------------------------
class PD_LPS_ADR_STRAIN : public PD_LPS_CRIT_STRAIN {
public:
  PD_LPS_ADR_STRAIN(PD_Particles &particles, double c, double stretchCrit,
                    double shearCrit, bool planeStress = false,
                    bool analyticalM = false);

  virtual void evaluateStatic(int id, int i);
  virtual void evaluateStepTwo(int id, int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_ADR_STRAIN_H
