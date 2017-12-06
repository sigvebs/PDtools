#ifndef PD_LPS_H
#define PD_LPS_H

#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{

//------------------------------------------------------------------------------
class PD_LPS : public Force
{
protected:
    bool m_planeStress;
//    int m_dim = 3;
    int m_iMicromodulus;
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
    int m_indexBrokenNow;
    int m_analyticalM;

    double m_k;
    double m_mu;
//    double m_delta;
//    double m_nu;
    double m_alpha;
    double m_c;
    double m_t;

    // Pointers
    double * m_mass;
    double * m_theta;
    double * m_volume;
    double * m_x;
    double * m_y;
    double * m_z;
    double * m_Fx;
    double * m_Fy;
    double * m_Fz;
    double * m_theta_new;

    int m_indexStress[6];
    vector<double*> m_stress;

    //       double weightFunction(const double dr0) const {return 1.;}
//           double weightFunction(const double dr0) const {return 1.0/dr0;}
        double weightFunction(const double dr0) const {return m_delta/dr0;}
//                   double weightFunction(const double dr0) const {return m_delta*m_delta/(dr0*dr0);}
//           double weightFunction(const double dr0) const {return exp(-3.*dr0/m_delta);}
    //       double weightFunction(const double dr0) const {return 1./(dr0*dr0);}
//    double weightFunction(const double dr0) const {return 1.*(1. - 1.01*dr0/m_delta);}
//           double weightFunction(const double dr0) const {return 0.5*(1. - 1.01*dr0/m_delta);}
public:
    PD_LPS(PD_Particles &particles, bool planeStress=false, bool analyticalM=false);

    virtual void
    calculateForces(const int id, const int i);

    virtual double
    calculatePotentialEnergyDensity(const int id_i, const int i);

    double
    computeDilation(const int id_i, const int i);

    virtual void
    calculatePotentialEnergy(const int id_i, const int i,
                                          int indexPotential);
    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);

    virtual void
    updateState(int id, int i);

    virtual double
    calculateStableMass(const int id_a, const int a,
                                     double dt);
    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    void
    calculateWeightedVolume();

    void
    updateWeightedVolume(int id_i, int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPS_H
