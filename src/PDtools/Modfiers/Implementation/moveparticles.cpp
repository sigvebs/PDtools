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
}
//------------------------------------------------------------------------------
MoveParticles::~MoveParticles()
{

}
//------------------------------------------------------------------------------
void MoveParticles::evaluateStepOne()
{
    arma::mat & r = m_particles->r();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;
        r(m_boundaryOrientation, col_i) += m_velAmplitude*m_dt;
//        for(int d=0;d<3; d++)
//        {
//            v(d, col_i) = 0.0;
//            F(d, col_i) = 0.0;
//            Fold(d, col_i) = 0.0;
//        }
    }
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
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double pos = r(m_boundaryOrientation, col_i);
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
        v(m_boundaryOrientation, col_i) = 0.0;
        F(m_boundaryOrientation, col_i) = 0.0;
        Fold(m_boundaryOrientation, col_i) = 0.0;

//        for(int d=0;d<3; d++)
//        {
//            v(d, col_i) = 0.0;
//            F(d, col_i) = 0.0;
//            Fold(d, col_i) = 0.0;
//        }
    }
}
//------------------------------------------------------------------------------
}
