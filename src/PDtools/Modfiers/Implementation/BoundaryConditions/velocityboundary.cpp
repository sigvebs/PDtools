#include "velocityboundary.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
VelocityBoundary::VelocityBoundary(double velAmplitude,
                                   double velOrientation,
                                   pair<double, double> boundary,
                                   int boundaryOrientation,
                                   double dt,
                                   int steps, int isStatic)
{
    m_velAmplitude = velAmplitude;
    m_velOritentation = velOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_dv = velAmplitude / steps;
    m_v = 0;
    m_dt = dt;
    m_isStatic = isStatic;

    if(m_velOritentation == 0)
    {
        m_otherAxis = {1, 2};
    }
    else if(m_velOritentation == 1)
    {
        m_otherAxis = {0, 2};
    }
    else if(m_velOritentation == 2)
    {
        m_otherAxis = {1, 2};
    }

    m_usingUnbreakableBorder = true;
}
//------------------------------------------------------------------------------
VelocityBoundary::~VelocityBoundary()
{

}
//------------------------------------------------------------------------------
void VelocityBoundary::registerParticleParameters()
{
    m_particles->registerParameter("unbreakable");
    m_ghostParameters = {"unbreakable"};
}
//------------------------------------------------------------------------------
void VelocityBoundary::evaluateStepOne()
{
    if(fabs(m_v) < fabs(m_velAmplitude))
    {
        m_v += m_dv;
    }
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & v = m_particles->v();
    arma::mat & r = m_particles->r();
    arma::mat & F = m_particles->F();
    arma::imat & isStatic = m_particles->isStatic();

    const double v_dt = m_v * m_dt;

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);

        v(i, m_velOritentation) = m_v;
        if(isStatic(i))
        {
            r(i, m_velOritentation) += v_dt;
        }
    }
}
//------------------------------------------------------------------------------
void VelocityBoundary::evaluateStepTwo()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::imat & isStatic = m_particles->isStatic();

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);

        v(i, m_velOritentation) = m_v;
        F(i, m_velOritentation) = 0;

        if(isStatic(i))
        {
            for(int d:m_otherAxis)
            {
                F(i, d) = 0;
            }
        }
    }
}
//------------------------------------------------------------------------------
void VelocityBoundary::initialize()
{
    // Selecting particles
    const arma::mat & r = m_particles->r();
    const ivec &colToId = m_particles->colToId();
    arma::mat & data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();
    const int unbreakablePos = m_particles->registerParameter("unbreakable");

    const double uRadius = 0.5*(m_boundary.second - m_boundary.first);
    int n = 0;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        const double pos = r(i, m_boundaryOrientation);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            const int id = colToId(i);

            if(m_isStatic)
                isStatic(i) = 1;
#ifdef USE_OPENMP
#pragma omp critical
#endif
            m_localParticleIds.push_back(id);
            data(i, unbreakablePos) = 1;
        }

        if(m_usingUnbreakableBorder)
        {
            if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
            {
#ifdef USE_OPENMP
#pragma omp critical
#endif
                data(i, unbreakablePos) = 1;
                n++;
            }
        }
    }
}
//------------------------------------------------------------------------------
}
