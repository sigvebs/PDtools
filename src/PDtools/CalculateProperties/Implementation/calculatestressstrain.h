#ifndef CALCULATESTRESSSTRAIN_H
#define CALCULATESTRESSSTRAIN_H

#include "PDtools/CalculateProperties/calculateproperty.h"
#include <vector>

namespace PDtools
{
class Force;
//------------------------------------------------------------------------------
class CalculateStressStrain : public CalculateProperty
{
public:
    CalculateStressStrain(vector<Force *> &forces,
                          double E, double nu, double delta,
                          bool planeStress=false,
                          bool greenLagrange=false);


    virtual void
    initialize();

    virtual void
    clean();

    virtual void
    update();

    void
    computeK(int id, int i);
private:
    vector<Force *> m_forces;
    int m_indexStress[6];
    int m_indexStrain[6];
    int m_indexK[6];
    int m_nStressStrainElements;

    int m_indexConnected;
    int m_indexVolume;
    int m_indexVolumeScaling;
    int m_indexDr0;
    int m_indexBrokenNow;
    bool m_smallStrain = false;
    bool m_greenStrain = false;
    bool m_planeStress = false;

    double m_E;
    double m_nu;
    double m_delta;
};
//------------------------------------------------------------------------------
}

#endif // CALCULATESTRESSSTRAIN_H
