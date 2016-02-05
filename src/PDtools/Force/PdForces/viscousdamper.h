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

    virtual void
    calculateForces(const int id, const int i);

protected:
    double m_c;
    double m_dt;
    arma::mat & m_v;
    arma::imat & m_isStatic;
    int m_indexDr0;
    int m_indexStretch_prev;
    int m_indexStretch;
    bool m_local;
};
//------------------------------------------------------------------------------
}
#endif // VISCOUSDAMPER_H
