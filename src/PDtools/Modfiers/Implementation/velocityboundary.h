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
                     int steps, bool isStatic);
    ~VelocityBoundary();

    virtual void evaluateStepOne();
    virtual void evaluateStepTwo();
    virtual void initialize();

private:
    vector<pair<int, int>> m_boundaryParticles;
    double m_velAmplitude;
    int m_velOritentation;
    pair<double, double> m_boundary;
    int m_boundaryOrientation;
    vector<int> m_otherAxis;
    double m_dv;
    double m_v;
    double m_dt;
    bool m_isStatic;
};
//------------------------------------------------------------------------------
}
#endif // VELOCITYBOUNDARY_H
