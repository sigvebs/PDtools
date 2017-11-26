#ifndef CALCULATEDAMAGE_H
#define CALCULATEDAMAGE_H

#include "PDtools/CalculateProperties/calculateproperty.h"

namespace PDtools
{
//------------------------------------------------------------------------------

class CalculateDamage : public CalculateProperty
{
public:
    CalculateDamage(double delta);

    virtual void initialize();

    virtual void update();

    double weightFunction(const double dr0) const {return m_delta/dr0;}

private:
    double m_delta;
    int m_iDamage;
    int m_iInitialWeight;
    int m_indexConnected;
    int m_iVolume;
    int m_iVolumeScaling;
    int m_iDr0;
};
//------------------------------------------------------------------------------
}

#endif // CALCULATEDAMAGE_H
