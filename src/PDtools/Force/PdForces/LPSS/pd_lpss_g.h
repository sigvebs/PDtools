#ifndef PD_LPSS_G_H
#define PD_LPSS_G_H

//#include <unordered_map>
#include "PDtools/Force/force.h"
#include "PDtools/Force/PdForces/LPSS/pd_lpss.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class PD_LPSS_G : public PD_LPSS
{
protected:
    double m_G0;
    double m_e_max;
    int m_iUnbreakable;
public:
    PD_LPSS_G(PD_Particles &particles, double G0, bool planeStress=false);

    virtual void
    calculateForces(const int id_i, const int i);

    virtual void
    evaluateStepTwo(const int id_i, const int i);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);
};
//------------------------------------------------------------------------------
}

#endif // PD_LPSS_G_H
