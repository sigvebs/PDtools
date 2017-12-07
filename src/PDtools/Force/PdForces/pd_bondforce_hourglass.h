#ifndef PD_BONDFORCE_HOURGLASS_H
#define PD_BONDFORCE_HOURGLASS_H

#include "PDtools/Force/PdForces/pd_bondforce.h"

namespace PDtools {

//------------------------------------------------------------------------------
class PD_bondForce_hourglass : public PD_bondForce {
public:
  PD_bondForce_hourglass(PD_Particles &particles);

  virtual void initialize(double E, double nu, double delta, int dim, double h,
                          double lc);

  virtual void calculateForces(const int id_i, const int i);

protected:
  int m_iF[6];
  int n_deformationGradient;
  int m_iC_hg;
};
//------------------------------------------------------------------------------
}
#endif // PD_BONDFORCE_HOURGLASS_H
