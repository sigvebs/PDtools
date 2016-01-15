#ifndef PD_BONDFORCE_H
#define PD_BONDFORCE_H

#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{

//------------------------------------------------------------------------------
class PD_bondForce : public Force
{
private:
    int m_indexMicromodulus;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexVolumeScaling;
    int m_indexStretch;
    int m_indexForceScaling;
    int m_indexWeightfunction;
    int m_indexConnected;
    int m_indexMyPdPosition;

    enum PD_bondForceErrorMessages
    {
        MicrmodulusNotSet
    };

public:
    PD_bondForce(PD_Particles &particles);

    ~PD_bondForce();

    virtual void
    calculateForces(const int id_i, const int i);

    virtual void
    calculateLinearForces(const int id_i, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

    virtual void
    calculatePotentialEnergy(const int id_i, const int i,
                                          int indexPotential);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);

    virtual double
    calculateStableMass(const int id_a, const int a, double dt);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);
};
//------------------------------------------------------------------------------
}
#endif // PD_BONDFORCE_H
