#ifndef ADR_H
#define ADR_H

#include <armadillo>
#include "solver.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class ADR: public Solver
{
protected:
    double m_globalError;
    double m_c;
    double m_epsilon = 0.2;
//    double m_errorThreshold = 0.001;

public:
    ADR();
    ~ADR();

    virtual void solve();

    virtual void stepForward(int i);
    void iterate();

protected:
    virtual void checkInitialization();
    virtual void initialize();
    void calculateStableMass();
    virtual void integrateStepOne();
    virtual void integrateStepTwo();
    void staticModifiers();
};
//------------------------------------------------------------------------------
}

#endif // ADR_H
