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
    double m_du_u;
    double m_c;
    double m_epsilon = 0.2;
    int m_maxSteps = 2000;
    int m_maxStepsFracture = 1000;
    int m_counter = 0;

public:
    ADR();

    virtual void
    solve();

    virtual void
    stepForward(int i);
    virtual void
    iterate(int maxNumberOfSteps=2500);

    void
    maxSteps(double _maxSteps){
        m_maxSteps = _maxSteps;
    }
    void
    maxStepsFracture(double _maxStepsFracture){
        m_maxStepsFracture = _maxStepsFracture;
    }

protected:
    virtual void
    checkInitialization();
    virtual void
    initialize();
    virtual void
    integrateStepOne();
    virtual void
    integrateStepTwo();
    void
    staticModifiers();
    void
    calculateStableMass();
    virtual void
    updateGridAndCommunication();
};
//------------------------------------------------------------------------------
}

#endif // ADR_H
