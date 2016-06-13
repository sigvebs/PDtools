#include "strainboundary.h"
#include "PDtools/Particles/pd_particles.h"

//#include <mpi.h>

namespace PDtools
{
//------------------------------------------------------------------------------
StrainBoundary::StrainBoundary(double strainrate,
                               double strainOrientation,
                               double nu,
                               pair<double, double> boundary,
                               int boundaryOrientation,
                               double dt,
                               bool isStatic)
{
    m_strainrate = strainrate;
    m_strainOritentation = strainOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_nu = nu;
    m_dt = dt;
    m_isStatic = isStatic;
    m_time = dt;
    m_usingUnbreakableBorder = true;
}
//------------------------------------------------------------------------------
void StrainBoundary::registerParticleParameters()
{
    m_particles->registerParameter("unbreakable");
    m_ghostParameters = {"unbreakable"};
}
//------------------------------------------------------------------------------
void StrainBoundary::evaluateStepOne()
{
    mat & r = m_particles->r();
    const mat & r0 = m_particles->r0();

    const double c1 = m_dt*m_strainrate;
    const double c2 = -m_nu*m_dt*m_strainrate;
    const int nParticles = m_particles->nParticles() + m_particles->nGhostParticles();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<nParticles; i++)
    {
        r(i, m_boundaryOrientation) += c1*r0(i, m_boundaryOrientation);

        for(int d:m_otherDirections)
        {
            r(i, d) += c2*r0(i, d);
        }
    }
}
//------------------------------------------------------------------------------
void StrainBoundary::initialize()
{
    m_particles->setNeedGhostR0(true);

    if(m_dim == 2)
    {
        switch (m_boundaryOrientation) {
        case 0:
            m_otherDirections = {1};
            break;
        case 1:
            m_otherDirections = {0};
            break;
        default:
            cerr << "ERROR for Modifier 'strain boundary:' strain in direction "
                 << m_boundaryOrientation
                 << " is not valid";
            exit(EXIT_FAILURE);
            break;
        }
    }

    if(m_dim == 3)
    {
        switch (m_boundaryOrientation) {
        case 0:
            m_otherDirections = {1, 2};
            break;
        case 1:
            m_otherDirections = {0, 2};
            break;
        case 2:
            m_otherDirections = {0, 1};
            break;
        default:
            cerr << "ERROR for Modifier 'strain boundary:' strain in direction "
                 << m_boundaryOrientation
                 << " is not valid";
            exit(EXIT_FAILURE);
            break;
        }
    }

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
void StrainBoundary::staticEvaluation()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & r = m_particles->r();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);

//        v(i, m_boundaryOrientation) = 0.0;
//        F(i, m_boundaryOrientation) = 0.0;
//        Fold(i, m_boundaryOrientation) = 0.0;

        for(int d=0;d<m_dim; d++)
        {
            v(i, d) = 0.0;
            F(i, d) = 0.0;
            Fold(i, d) = 0.0;
        }
    }
}
//------------------------------------------------------------------------------
}
