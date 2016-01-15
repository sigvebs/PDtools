#ifndef VISCOUSDAMPER_H
#define VISCOUSDAMPER_H

#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class ViscousDamper : public Force
{
public:
    ViscousDamper(PD_Particles &particles, const double dampeningCoeff);
    ~ViscousDamper();

    virtual void
    calculateForces(const int id, const int i);

protected:
    double m_c;
    arma::mat & m_v;
    arma::imat & m_isStatic;
};
//------------------------------------------------------------------------------
}
#endif // VISCOUSDAMPER_H
