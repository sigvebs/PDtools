#ifndef BOUNDARYFORCE_H
#define BOUNDARYFORCE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class boundaryForce : public Modifier
{
public:
    boundaryForce(double appliedForce, double forceOrientation,
                     pair<double, double> boundary, int boundaryOrientation, int steps, double delta);
    ~boundaryForce();

    virtual void evaluateStepOne();
    virtual void evaluateStepTwo();
    virtual void initialize();
    virtual void staticEvaluation();

private:
    vector<pair<int, int>> m_boundaryParticles;
    double m_forceDensity;
    int m_forceOritentation;
    pair<double, double> m_boundary;
    int m_boundaryOrientation;
    vector<int> m_otherAxis;
    int m_indexVolume ;
    int m_steps;
    int m_indexRadius;
    double m_incrementalForce;
    int m_inceremental;
    double m_delta;
};
//------------------------------------------------------------------------------
}
#endif // BOUNDARYFORCE_H
