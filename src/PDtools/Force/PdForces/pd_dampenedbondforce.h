#ifndef PD_DAMPENEDBONDFORCE_H
#define PD_DAMPENEDBONDFORCE_H

#include "PDtools/Force/PdForces/pd_bondforce.h"

namespace PDtools
{

//------------------------------------------------------------------------------
class PD_dampenedBondForce : public PD_bondForce
{
public:
    PD_dampenedBondForce(PD_Particles &particles, double dt, double c);

    virtual void
    calculateForces(const int id_i, const int i);

protected:
    double m_dt;
    double m_c;
};
//------------------------------------------------------------------------------
}
#endif // PD_DAMPENEDBONDFORCE_H
