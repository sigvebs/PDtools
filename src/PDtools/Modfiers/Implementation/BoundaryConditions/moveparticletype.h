#ifndef MOVEPARTICLETYPE_H
#define MOVEPARTICLETYPE_H


#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MoveParticleGroup : public Modifier
{
public:
    MoveParticleGroup(int type, double velAmplitude, vec velocityDirection,
                      double dt, bool isStatic);

    virtual void
    registerParticleParameters();

    virtual void
    evaluateStepOne();

    virtual void
    staticEvaluation();

    virtual void
    initialize();

private:
    int m_groupId;
    double m_velAmplitude;
    double m_dt;
    bool m_isStatic;
    double m_time;
    vec m_velocityDirection;
    vec m_vd = {0,0,0};
};
//------------------------------------------------------------------------------
}
#endif // MOVEPARTICLETYPE_H
