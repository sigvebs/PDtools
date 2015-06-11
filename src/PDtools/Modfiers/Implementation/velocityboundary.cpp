#include "velocityboundary.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
VelocityBoundary::VelocityBoundary(double velAmplitude,
                                   double velOrientation,
                                   pair<double, double> boundary,
                                   int boundaryOrientation,
                                   int steps)
{
    m_velAmplitude = velAmplitude;
    m_velOritentation = velOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_dv = velAmplitude / steps;
    m_v = 0;
}
//------------------------------------------------------------------------------
VelocityBoundary::~VelocityBoundary()
{

}
//------------------------------------------------------------------------------
void VelocityBoundary::evaluateStepTwo()
{
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;
        v(m_boundaryOrientation, col_i) = m_v;
        F(m_boundaryOrientation, col_i) = 0.0;
    }

    if(fabs(m_v) < fabs(m_velAmplitude))
    {
        m_v += m_dv;
    }
}
//------------------------------------------------------------------------------
void VelocityBoundary::initialize()
{
    // Selecting particles
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();
    int unbreakablePos;

    if(m_particles->hasParameter("unbreakable"))
    {
        unbreakablePos = m_particles->getParamId("unbreakable");
    }
    else
    {
        unbreakablePos = m_particles->registerParameter("unbreakable");
    }

    double uRadius = 0.5*(m_boundary.second - m_boundary.first);
    int n = 0;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double pos = r(m_boundaryOrientation, col_i);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            pair<int, int> pId(i, i);
#pragma omp critical
            m_boundaryParticles.push_back(pId);
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
        {
//#pragma omp critical
            data(i, unbreakablePos) = 1;
            n++;
        }
    }
}
//------------------------------------------------------------------------------
}
