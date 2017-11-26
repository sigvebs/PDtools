#ifndef PD_PMB_H
#define PD_PMB_H

#include <unordered_map>
#include "PDtools/Force/force.h"

//------------------------------------------------------------------------------
namespace PDtools
{

class PD_PMB : public Force
{
public:
    int m_dim = 3;
    int m_indexMicromodulus;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexVolumeScaling;
    int m_indexStretch;
    int m_indexForceScaling;
    int m_indexConnected;
    double m_lc;
    double m_delta;
    int m_indexS0;
    int m_indexS_new;
    int m_indexS00;
    int m_indexUnbreakable;

    double **f;
    double **x;
    double **r0;

    double m_alpha;
    int m_indexStress[6];

    enum PD_bondForceErrorMessages {
        MicrmodulusNotSet
    };

public:
    PD_PMB(PD_Particles &particles, double lc, double delta, double alpha);
    ~PD_PMB();

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

    virtual double
    calculateStableMass(const int id_a, const int a,
                                     double dt);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);
};
//------------------------------------------------------------------------------
}

#endif // PD_PMB_H
