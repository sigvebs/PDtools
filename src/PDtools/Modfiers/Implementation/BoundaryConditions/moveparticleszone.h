#ifndef MOVEPARTICLESZONE_H
#define MOVEPARTICLESZONE_H


#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MoveParticlesZone : public Modifier
{
public:
    MoveParticlesZone(double velAmplitude, vec velocityDirection,
                      vector<double> initialArea,
                      double dt, bool isStatic, double delta);

    virtual void
    registerParticleParameters();

    virtual void
    evaluateStepOne();

    virtual void
    staticEvaluation();

    virtual void
    initialize();

private:
    double m_delta;
    double m_velAmplitude;
    vector<double> m_initialArea;
    double m_dt;
    bool m_isStatic;
    double m_time;
    bool m_usingUnbreakableBorder;
    vec m_velocityDirection;
    vec m_vd_abs = {0,0,0};
};
//------------------------------------------------------------------------------
}
#endif // MOVEPARTICLESZONE_H
