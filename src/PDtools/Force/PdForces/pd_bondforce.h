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
    int m_indexCompute;

    enum PD_bondForceErrorMessages
    {
        MicrmodulusNotSet
    };

public:
    PD_bondForce(PD_Particles &particles);

    ~PD_bondForce();

    virtual void
    calculateForces(const std::pair<int, int> & idCol);

    virtual void
    calculateLinearForces(const std::pair<int, int> & idCol);

    virtual double
    calculatePotentialEnergyDensity(const std::pair<int, int> & idCol);

    virtual void
    calculatePotentialEnergy(const std::pair<int, int> & idCol,
                                          int indexPotential);

    virtual void
    calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);

    virtual double
    calculateStableMass(const std::pair<int, int> & idCol,
                                     double dt);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);
};
//------------------------------------------------------------------------------
}
#endif // PD_BONDFORCE_H
