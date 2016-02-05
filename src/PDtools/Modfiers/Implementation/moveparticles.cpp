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
    m_time = dt;
    m_usingUnbreakableBorder = true;
}
//------------------------------------------------------------------------------
MoveParticles::~MoveParticles()
{

}
//------------------------------------------------------------------------------
void MoveParticles::registerParticleParameters()
{
    m_particles->registerParameter("unbreakable");
    m_ghostParameters = {"unbreakable"};
}
//------------------------------------------------------------------------------
void MoveParticles::evaluateStepOne()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    mat & r = m_particles->r();
//    const mat & r0 = m_particles->r0();
    const double dr = m_time*m_velAmplitude;
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i =  idToCol.at(id);
//        for(int d=0;d<DIM; d++)
//        {
//            r(i, d) = r0(i, d);
//        }
//        r(i, m_boundaryOrientation) = r0(i, m_boundaryOrientation) + dr;
        r(i, m_boundaryOrientation) += dr;

        for(int d=0;d<DIM; d++)
        {
            v(i, d) = 0.0;
            F(i, d) = 0.0;
            Fold(i, d) = 0.0;
        }
    }
//    cout << dr << " " << m_time << " " << m_velAmplitude << " " << m_dt <<  endl;
//    m_time += m_dt;
}
//------------------------------------------------------------------------------
void MoveParticles::initialize()
{
    // Selecting particles
    const ivec &colToId = m_particles->colToId();
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();
    const int unbreakablePos = m_particles->registerParameter("unbreakable");

    double uRadius = 0.5*(m_boundary.second - m_boundary.first);
    int n = 0;

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        double pos = r(i, m_boundaryOrientation);
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
void MoveParticles::staticEvaluation()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);
        for(int d=0;d<DIM; d++)
        {
            v(i, d) = 0.0;
            F(i, d) = 0.0;
            Fold(i, d) = 0.0;
        }
    }
}
//------------------------------------------------------------------------------
}
