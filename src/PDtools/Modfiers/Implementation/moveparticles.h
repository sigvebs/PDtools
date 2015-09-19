#ifndef MOVEPARTICLES_H
#define MOVEPARTICLES_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MoveParticles: public Modifier
{
public:
    MoveParticles(double velAmplitude, double velOrientation,
                   pair<double, double> boundary, int boundaryOrientation,
                  double dt, bool isStatic);
    ~MoveParticles();

    virtual void evaluateStepOne();
    virtual void initialize();
    virtual void staticEvaluation();

private:
    vector<pair<int, int>> m_boundaryParticles;
    double m_velAmplitude;
    int m_velOritentation;
    pair<double, double> m_boundary;
    int m_boundaryOrientation;
    double m_dt;
    bool m_isStatic;
};
//------------------------------------------------------------------------------
}

#endif // MOVEPARTICLES_H
