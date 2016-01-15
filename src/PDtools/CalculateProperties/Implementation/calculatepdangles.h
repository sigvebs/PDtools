#ifndef CALCULATEPDANGLES_H
#define CALCULATEPDANGLES_H

#include "PDtools/CalculateProperties/calculateproperty.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class CalculatePdAngles : public CalculateProperty
{
public:
    CalculatePdAngles();
    ~CalculatePdAngles();

    virtual void
    initialize();

    virtual void
    update();
protected:
    int m_indexTheta;
    int m_indexTheta0;
    int m_indexConnected;
};
//------------------------------------------------------------------------------
}
#endif // CALCULATEPDANGLES_H
