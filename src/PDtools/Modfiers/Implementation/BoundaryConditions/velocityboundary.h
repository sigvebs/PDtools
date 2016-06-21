#ifndef VELOCITYBOUNDARY_H
#define VELOCITYBOUNDARY_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class VelocityBoundary: public Modifier
{
public:
    VelocityBoundary(double velAmplitude, double velOrientation,
                     pair<double, double> boundary, int boundaryOrientation, double dt,
                     int steps, int isStatic);
    ~VelocityBoundary();

    virtual void registerParticleParameters();
    virtual void evaluateStepOne();
    virtual void evaluateStepTwo();
    virtual void initialize();

private:
    double m_velAmplitude;
    int m_velOritentation;
    pair<double, double> m_boundary;
    int m_boundaryOrientation;
    vector<int> m_otherAxis;
    double m_dv;
    double m_v;
    double m_dt;
    bool m_isStatic;
    bool m_usingUnbreakableBorder = false;
};
//------------------------------------------------------------------------------
}
#endif // VELOCITYBOUNDARY_H
