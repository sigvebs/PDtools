#include "boundaryforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
boundaryForce::boundaryForce(double appliedForce,
                                   double forceOrientation,
                                   pair<double, double> boundary,
                                   int boundaryOrientation)
{
    m_forceDensity = appliedForce;
    m_forceOritentation = forceOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
}
//------------------------------------------------------------------------------
boundaryForce::~boundaryForce()
{

}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepOne()
{
    arma::mat & b = m_particles->b();
    const arma::mat & data = m_particles->data();
    arma::mat & F = m_particles->F();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        const int col_i = idCol.second;
        const double V_i = data(col_i, m_indexVolume);
        const double force = m_forceDensity;
//        b(col_i, m_forceOritentation) = force;
        F(col_i, m_forceOritentation) += force;
    }
}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepTwo()
{
    arma::mat & b = m_particles->b();
    b.zeros();
    evaluateStepOne();
}
//------------------------------------------------------------------------------
void boundaryForce::staticEvaluation()
{
    boundaryForce::evaluateStepTwo();
}
//------------------------------------------------------------------------------
void boundaryForce::initialize()
{
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();

    m_indexVolume = m_particles->getParamId("volume");
    double volume = 0;

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

    // Selecting particles
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        const int col_i = i;
        const double pos = r(col_i, m_boundaryOrientation);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            const pair<int, int> pId(i, i);
            m_boundaryParticles.push_back(pId);
            volume += data(i, m_indexVolume);
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
        {
            data(i, unbreakablePos) = 1;
        }
    }
    m_forceDensity /= volume;
}
//------------------------------------------------------------------------------
}
