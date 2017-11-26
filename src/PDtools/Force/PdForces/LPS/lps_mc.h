#ifndef LPS_MC_H
#define LPS_MC_H

#include<armadillo>
#include "PDtools/Force/PdForces/LPS/pd_lps.h"

using namespace arma;

namespace PDtools
{
//------------------------------------------------------------------------------
class LPS_mc : public PD_LPS
{
public:
    LPS_mc(PD_Particles &particles, double c, double phi, double C, double T, bool planeStress=false, bool analyticalM=false);

    virtual void
    calculateForces(const int id, const int i);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    void
    computeK(int id, int i);

    void
    updateWeightedMass(int id, int i);

    virtual void
    evaluateStepTwo(int id, int i);
protected:
    void computeStress(const int id, const int i, const int nConnected);

    // Material dampening
    double m_dampCoeff;

    // MC criterion
    double m_phi;
    double m_d;
    double m_S0;
    double m_T;
    double m_cos_theta;
    double m_sin_theta;
    double m_tan_theta;
    double m_ks;
    double m_C0;

    // Stress strain related
    int m_indexStress[6];
    int m_indexStrain[6];
    int m_indexK[6];
    int m_nStressStrainElements;

    int m_indexConnected;
    int m_indexUnbreakable;
    int m_indexVolume;
    int m_indexVolumeScaling;
    int m_indexDr0;
    int m_indexBrokenNow;
    bool m_smallStrain = false;
    bool m_greenStrain = false;
    bool m_planeStress = false;

    mat _F;
    mat m_strain;
    mat m_K;
    mat m_P;

    double m_lambda;
    double m_mu;
};
//------------------------------------------------------------------------------
}

#endif // LPS_MC_H
