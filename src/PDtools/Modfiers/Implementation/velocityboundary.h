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
                     pair<double, double> boundary, int boundaryOrientation,
                     int steps=1);
    ~VelocityBoundary();

    void evaluateStepTwo();
    virtual void initialize();

private:
    vector<pair<int, int>> m_boundaryParticles;
    double m_velAmplitude;
    int m_velOritentation;
    pair<double, double> m_boundary;
    int m_boundaryOrientation;
    double m_dv;
    double m_v;
};
//------------------------------------------------------------------------------
}
#endif // VELOCITYBOUNDARY_H
