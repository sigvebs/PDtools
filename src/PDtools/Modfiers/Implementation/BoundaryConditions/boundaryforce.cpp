#include "boundaryforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
boundaryForce::boundaryForce(double appliedForce,
                                   double forceOrientation,
                                   pair<double, double> boundary,
                                   int boundaryOrientation,
                                   int steps,
                             double delta)
{
    m_forceDensity = appliedForce;
    m_forceOritentation = forceOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_inceremental = false;
    m_incrementalForce = 0;
    m_delta = delta;
    if(!m_inceremental)
    {
        m_incrementalForce = m_forceDensity;
    }
}
//------------------------------------------------------------------------------
boundaryForce::~boundaryForce()
{

}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepOne()
{
    arma::mat & F = m_particles->F();
    arma::mat & data = m_particles->data();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        const int i = idCol.second;
        const double radius = data(i, m_indexRadius);
        F(i, m_forceOritentation) += m_incrementalForce/radius;
    }
}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepTwo()
{
    arma::mat & b = m_particles->b();
    b.zeros();
    if(m_inceremental)
    {
        m_incrementalForce += m_forceDensity;
    }
    evaluateStepOne();
}
//------------------------------------------------------------------------------
void boundaryForce::staticEvaluation()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);

        v(i, m_forceOritentation) = 0.0;
        F(i, m_forceOritentation) = 0.0;
        Fold(i, m_forceOritentation) = 0.0;

//        for(int d=0;d<M_DIM; d++)
//        {
//            v(i, d) = 0.0;
//            F(i, d) = 0.0;
//            Fold(i, d) = 0.0;
//        }
    }
//    boundaryForce::evaluateStepTwo();
}
//------------------------------------------------------------------------------
void boundaryForce::initialize()
{
    m_indexRadius = m_particles->registerParameter("radius");
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

//    double uRadius = 0.5*(m_boundary.second - m_boundary.first);
    double uRadius = 6*m_delta;

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
            data(i, unbreakablePos) = 1;
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
        {
            data(i, unbreakablePos) = 1;
        }
    }
//    m_forceDensity /= volume;
}
//------------------------------------------------------------------------------
}
