#ifndef PD_LPSDAMPENEDCONTACT_P_H
#define PD_LPSDAMPENEDCONTACT_P_H

#include "PDtools/Force/PdForces/LPS_porosity/pd_lps_p.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class PD_lpsDampenedContact_porosity : public PD_LPS_POROSITY
{
public:
    PD_lpsDampenedContact_porosity(PD_Particles &particles, double m, double b, double c, bool planeStress);

    virtual void
    calculateForces(const int id, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

protected:
    double m_dampCoeff;
};
//------------------------------------------------------------------------------
}
#endif // PD_LPSDAMPENEDCONTACT_P_H
