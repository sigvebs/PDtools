#ifndef DEMFORCE_H
#define DEMFORCE_H


#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class DemForce : public Force
{
public:
    DemForce(PD_Particles &particles, double T, double C);

    virtual void
    calculateForces(const int id_i, const int i);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);
private:
    int m_indexMicromodulus;
    int m_indexShearModulus;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexVolumeScaling;
    int m_indexStretch;
    int m_indexConnected;
    int m_indexRadius;
    int m_indexUnbreakable;
    int m_indexFs[M_DIM];

    double m_Kn;
    double m_Ks;
    double m_T;
    double m_C;
};
//------------------------------------------------------------------------------
}
#endif // DEMFORCE_H
