#ifndef PD_LPSDAMPENEDCONTACT_H
#define PD_LPSDAMPENEDCONTACT_H

#include "PDtools/Force/PdForces/pd_lps.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class PD_lpsDampenedContact : public PD_LPS
{
public:
    PD_lpsDampenedContact(PD_Particles &particles, double c, bool planeStress);

    virtual void
    calculateForces(const int id, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

protected:
    double m_dampCoeff;
};
//------------------------------------------------------------------------------
}
#endif // PD_LPSDAMPENEDCONTACT_H
