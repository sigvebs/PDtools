#include "moveparticles.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
MoveParticles::MoveParticles(double velAmplitude,
                             double velOrientation,
                             pair<double, double> boundary,
                             int boundaryOrientation,
                             double dt,
                             bool isStatic)
{
    m_velAmplitude = velAmplitude;
    m_velOritentation = velOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_dt = dt;
    m_isStatic = isStatic;
    m_time = 0;
}
//------------------------------------------------------------------------------
MoveParticles::~MoveParticles()
{

}
//------------------------------------------------------------------------------
void MoveParticles::evaluateStepOne()
{
    mat & r = m_particles->r();
    const mat & r0 = m_particles->r0();
    const double dr = m_time*m_velAmplitude;
//    mat & v = m_particles->v();
//    mat & F = m_particles->F();
//    mat & Fold = m_particles->Fold();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;
        for(int d=0;d<3; d++)
        {
            r(col_i, d) = r0(col_i, d);
        }
        r(col_i, m_boundaryOrientation) = r0(col_i, m_boundaryOrientation) + dr;
    }

    m_time += m_dt;
}
//------------------------------------------------------------------------------
void MoveParticles::initialize()
{
    // Selecting particles
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();
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
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double pos = r(col_i, m_boundaryOrientation);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            pair<int, int> pId(i, i);
            if(m_isStatic)
                isStatic(pId.second) = 1;
#ifdef USE_OPENMP
#pragma omp critical
#endif
            m_boundaryParticles.push_back(pId);
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
        {
            data(i, unbreakablePos) = 1;
            n++;
        }
    }
}
//------------------------------------------------------------------------------
void MoveParticles::staticEvaluation()
{
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;
//        v(col_i, m_boundaryOrientation) = 0.0;
//        F(col_i, m_boundaryOrientation) = 0.0;
//        Fold(col_i, m_boundaryOrientation) = 0.0;
        for(int d=0;d<3; d++)
        {
            v(col_i, d) = 0.0;
            F(col_i, d) = 0.0;
            Fold(col_i, d) = 0.0;
        }
    }
}
//------------------------------------------------------------------------------
}
