#ifndef CALCULATEDAMAGE_H
#define CALCULATEDAMAGE_H

#include "PDtools/CalculateProperties/calculateproperty.h"

namespace PDtools
{
//------------------------------------------------------------------------------

class CalculateDamage : public CalculateProperty
{
public:
    CalculateDamage();

    virtual void
    initialize();

    virtual void
    update();

private:
    int m_indexDamage;
    int m_indexMaxPdConnections;
    int m_indexConnected;
};
//------------------------------------------------------------------------------
}

#endif // CALCULATEDAMAGE_H
