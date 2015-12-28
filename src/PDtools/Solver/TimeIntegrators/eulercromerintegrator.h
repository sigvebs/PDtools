#ifndef EULERCROMERINTEGRATOR_H
#define EULERCROMERINTEGRATOR_H


#include "PDtools/Solver/timeintegrator.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class EulerCromerIntegrator : public TimeIntegrator
{
public:
    EulerCromerIntegrator();
    ~EulerCromerIntegrator();

    virtual void integrateStepOne();
    virtual void integrateStepTwo();
};
//------------------------------------------------------------------------------
}
#endif // EULERCROMERINTEGRATOR_H
