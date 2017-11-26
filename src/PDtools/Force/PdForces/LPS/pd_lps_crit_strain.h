#ifndef PD_LPS_CRIT_STRAIN_H
#define PD_LPS_CRIT_STRAIN_H


#include<armadillo>
#include "PDtools/Force/PdForces/LPS/pd_lps.h"

using namespace arma;
namespace PDtools
{
//------------------------------------------------------------------------------
class PD_LPS_CRIT_STRAIN : public PD_LPS
{
public:
    PD_LPS_CRIT_STRAIN(PD_Particles &particles, double c, double stretchCrit, double shearCrit,
                       bool planeStress=false, bool analyticalM=false);


    virtual void
    evaluateStepTwo(int id_i, int i);
protected:
    double m_stretchCrit;
    double m_shearCrit;

    // Material dampening
    double m_dampCoeff;

    int m_iUnbreakable;
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_CRIT_STRAIN_H
