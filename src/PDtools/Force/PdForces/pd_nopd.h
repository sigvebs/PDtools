#ifndef PD_NOPD_H
#define PD_NOPD_H


//#include <unordered_map>
#include<armadillo>
#include "PDtools/Force/force.h"

using namespace arma;

namespace PDtools
{
//------------------------------------------------------------------------------
class PD_NOPD : public Force
{
public:
    PD_NOPD(PD_Particles &particles, double phi, double C, double T, bool planeStress=true);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    virtual void
    calculateForces(const int id, const int i);

    virtual void
    updateState(int id, int i);

    virtual void
    evaluateStepTwo(int id, int i);

    void
    computeK(int id, int i);
protected:

    // Material dampening
    double m_dampCoeff;

    // MC criterion
    double m_phi;
    double m_d;
    double m_C;
    double m_T;
    double m_cos_theta;
    double m_sin_theta;

    // Stress strain related
    int m_indexStress[6];
    int m_indexStrain[6];
    int m_indexK[6];
    int m_nStressStrainElements;

    int m_iConnected;
    int m_iUnbreakable;
    int m_iVolume;
    int m_iVolumeScaling;
    int m_iDr0;
    int m_iBrokenNow;
    bool m_smallStrain = false;
    bool m_greenStrain = false;
    bool m_planeStress = false;


    mat _F;
    mat m_strain;
    mat m_K_i;
    mat m_K_j;
    mat m_P_i;
    mat m_P_j;
    vec f_ij;
};
//------------------------------------------------------------------------------
}

#endif // PD_NOPD_H
