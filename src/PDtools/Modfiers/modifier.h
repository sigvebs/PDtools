#ifndef MODIFIER_H
#define MODIFIER_H

#include <armadillo>

using namespace std;

namespace PDtools
{
// Forward declerations
class PD_Particles;

//------------------------------------------------------------------------------
class Modifier
{
protected:
    PD_Particles *m_particles;
    bool m_state = false;

public:
    Modifier();
    virtual ~Modifier();

    virtual void evaluateStepOne(const pair<int, int> &id_col);
    virtual void evaluateStepTwo(const pair<int, int> &id_col);
    virtual void evaluateStepOne();
    virtual void evaluateStepTwo();
    virtual void staticEvaluation();
    virtual void initialize();
    void setParticles(PD_Particles &particles)
    {
        m_particles = &particles;
    }

    bool state();
};
//------------------------------------------------------------------------------
}
#endif // MODIFIER_H
