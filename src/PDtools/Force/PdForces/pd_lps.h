#ifndef PD_LPS_H
#define PD_LPS_H

#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{

//------------------------------------------------------------------------------
class PD_LPS : public Force
{
private:
    int m_dim = 3;
    int m_iTheta;
    int m_iThetaNew;
    int m_iVolume;
    int m_iDr0;
    int m_iVolumeScaling;
    int m_iStretch;
    int m_iForceScalingDilation;
    int m_iForceScalingBond;
    int m_iConnected;
    int m_iMass;
    int m_iCompute;

    arma::mat & m_r;
    arma::mat & m_r0;
    arma::mat & m_F;
    arma::mat & m_data;
    std::unordered_map<int, int> & m_pIds;

    double m_k;
    double m_mu;
    double m_delta;
    double m_nu;
public:
    PD_LPS(PD_Particles &particles);
    ~PD_LPS();
    virtual void calculateForces(const std::pair<int, int> & idCol);

    virtual double calculatePotentialEnergyDensity(const std::pair<int, int> & idCol);

    double computeDilation(const std::pair<int, int> & idCol);

    virtual void calculatePotentialEnergy(const std::pair<int, int> & idCol,
                                          int indexPotential);
    virtual void calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);

    virtual void updateState(const std::pair<int, int> & idCol);

    virtual double calculateStableMass(const std::pair<int, int> & idCol,
                                     double dt);
    virtual void initialize(double E, double nu, double delta, int dim, double h);

    void calculateMass();

    virtual void applySurfaceCorrection(double strain=0.01);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_H
