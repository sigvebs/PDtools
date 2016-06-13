#ifndef CALCULATESTRAIN_H
#define CALCULATESTRAIN_H

#include "PDtools/CalculateProperties/calculateproperty.h"
#include <vector>

namespace PDtools
{
//------------------------------------------------------------------------------
class CalculateStrain : public CalculateProperty
{
public:
    CalculateStrain(double delta, vector<pair<double, double> > &domain);

    virtual void
    initialize();

    virtual void
    clean();

    virtual void
    update();

    void
    calulateShapeFunction();

private:
    double m_delta;
    vector<pair<double,double>> &m_domain;
    int m_indexStrain[6];
    int m_indexShapeFunction[6];
    int m_nStrainElements;
    int m_indexConnected;
    int m_indexVolume;
    int m_indexVolumeScaling;
    int m_indexDr0;
    double m_w;
    double m_h = 1.;
};
//------------------------------------------------------------------------------
}

#endif // CALCULATESTRAIN_H
