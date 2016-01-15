#ifndef PD_OSP_H
#define PD_OSP_H

#include "PDtools/Force/force.h"

//------------------------------------------------------------------------------
namespace PDtools
{

class PD_OSP : public Force
{
private:
    int m_dim = 3;
    int m_indexA;
    int m_indexB;
    int m_indexD;
    int m_indexTheta;
    int m_indexThetaNew;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexVolumeScaling;
    int m_indexStretch;
    int m_indexForceScalingDilation;
    int m_indexForceScalingBond;
    int m_indexConnected;

    double m_delta;
public:
    PD_OSP(PD_Particles &particles);
    ~PD_OSP();

    virtual void
    calculateForces(const int id, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

    virtual void
    calculatePotentialEnergy(const int id_i, const int i,
                                          int indexPotential);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);

    virtual void
    updateState(int id, int i);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    virtual void
    applySurfaceCorrectionStep1(double mu, double nu, int dim, double strain);

    double
    calculateDilationTerm(const int id_i, const int i);

    double
    calculateBondPotential(const int id_i, const int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_OSP_H
