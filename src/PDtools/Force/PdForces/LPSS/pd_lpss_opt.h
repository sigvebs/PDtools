#ifndef PD_LPSS_OPT_H
#define PD_LPSS_OPT_H

#include "PDtools/Force/PdForces/LPSS/pd_lpss.h"

namespace PDtools {
//------------------------------------------------------------------------------
class PD_LPSS_opt : public PD_LPSS {
protected:
  int m_iRn[9]; // New Rotation matrix
  int n_nElements = 0;

public:
  PD_LPSS_opt(PD_Particles &particles, bool planeStress = false);

  virtual void calculateForces(const int id, const int i);
  virtual void evaluateStepOne();
  virtual void evaluateStepOne(const int id_i, const int i);
};
//------------------------------------------------------------------------------
}

#endif // PD_LPSS_OPT_H
