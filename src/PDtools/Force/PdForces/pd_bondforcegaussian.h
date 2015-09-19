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
    int m_indexExponent;
    int m_indexConnected;

    arma::mat & m_r;
    arma::mat & m_F;
    arma::mat & m_data;
    std::unordered_map<int, int> & m_pIds;
    double m_l;

    enum PD_bondForceErrorMessages
    {
        MicrmodulusNotSet
    };

public:
    PD_bondforceGaussian(PD_Particles &particles);
    ~PD_bondforceGaussian();
    virtual void calculateForces(const std::pair<int, int> & idCol);
    virtual double calculatePotentialEnergyDensity(const std::pair<int, int> & idCol);
    virtual void calculatePotentialEnergy(const std::pair<int, int> & idCol,
                                          int indexPotential);
    virtual double calculateBondEnergy(const std::pair<int, int> &idCol,
                                    std::pair<int, std::vector<double> > &con);
    virtual void calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);

    virtual double calculateStableMass(const std::pair<int, int> & idCol,
                                     double dt);

    virtual void initialize(double E, double nu, double delta, int dim, double h);
};
//------------------------------------------------------------------------------
}
#endif // PD_BONDFORCEGAUSSIAN_H
