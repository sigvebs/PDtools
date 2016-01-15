#include "calculatestress.h"

#include "Particles/pd_particles.h"
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculateStress::CalculateStress(vector<Force *> &forces):
    CalculateProperty("stress")
{
    for(Force* force: forces)
    {
        m_forces.push_back(force);
    }
}
//------------------------------------------------------------------------------

CalculateStress::~CalculateStress()
{

}
//------------------------------------------------------------------------------
void CalculateStress::initialize()
{
    switch(m_dim)
    {
    case 1:
        m_nStressElements = 1;
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        break;
    case 2:
        m_nStressElements = 3;
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        break;
    case 3:
        m_nStressElements = 6;
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStress[3] = m_particles->registerParameter("s_zz");
        m_indexStress[4] = m_particles->registerParameter("s_xz");
        m_indexStress[5] = m_particles->registerParameter("s_yz");
        break;
    }
}
//------------------------------------------------------------------------------
void CalculateStress::clean()
{
    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();// + m_particles->nGhostParticles();
    arma::mat & data = m_particles->data();

    // Zeroing out the stress
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        for(int s=0; s<m_nStressElements; s++)
        {
            data(i, m_indexStress[s]) = 0;
        }
    }
}
//------------------------------------------------------------------------------
void CalculateStress::update()
{
    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();
//    arma::mat & data = m_particles->data();

    // Updating single particle states
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        for(Force *force: m_forces)
        {
            force->calculateStress(id_i, i, m_indexStress);
        }
    }
}
//------------------------------------------------------------------------------
}
