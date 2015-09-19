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
void boundaryForce::evaluateStepTwo()
{
    arma::mat & F = m_particles->F();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;
        F(m_forceOritentation, col_i) += m_forceDensity;
    }
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

    const int indexVolume = m_particles->getParamId("volume");
    double volume = 0;

    // Selecting particles
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double pos = r(m_boundaryOrientation, col_i);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            pair<int, int> pId(i, i);
//#pragma omp critical
            m_boundaryParticles.push_back(pId);
            volume += data(i, indexVolume);
        }
    }
    m_forceDensity /= volume;
}
//------------------------------------------------------------------------------
}
