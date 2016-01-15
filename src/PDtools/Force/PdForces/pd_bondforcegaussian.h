#ifndef PD_BONDFORCEGAUSSIAN_H
#define PD_BONDFORCEGAUSSIAN_H

#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class PD_bondforceGaussian: public Force
{
private:
    int m_indexMicromodulus;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexVolumeScaling;
    int m_indexForceScaling;
    int m_indexStretch;
    int m_indexWeightFunction;
    int m_indexConnected;

    std::string m_weightType;

    double m_l;

    enum PD_bondForceErrorMessages
    {
        MicrmodulusNotSet
    };

public:
    PD_bondforceGaussian(PD_Particles &particles, std::string weightType="constant");
    ~PD_bondforceGaussian();

    virtual void
    calculateForces(const int id, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

    virtual void
    calculatePotentialEnergy(const int id_i, const int i,
                                          int indexPotential);

    virtual double
    calculateBondEnergy(const int id_i, const int i,
                                    std::pair<int, std::vector<double> > &con);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);

    virtual double
    calculateStableMass(const int id_a, const int a,
                                     double dt);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    void
    initializeConstant();

    void
    initializeLinear();

    void
    initializeGaussian();

    void
    initializeSigmoid();

};
//------------------------------------------------------------------------------
}
#endif // PD_BONDFORCEGAUSSIAN_H
