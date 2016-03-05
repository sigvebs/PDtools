#include "calculatedamage.h"

#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculateDamage::CalculateDamage():
    CalculateProperty("damage")
{

}
//------------------------------------------------------------------------------
void CalculateDamage::initialize()
{
    m_indexDamage = m_particles->registerParameter("damage");
    m_indexMaxPdConnections = m_particles->registerParameter("maxPdConnections");
    m_indexConnected = m_particles->getPdParamId("connected");

    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();
    mat & data = m_particles->data();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        data(i, m_indexMaxPdConnections) = m_particles->pdConnections(id_i).size();
    }
}
//------------------------------------------------------------------------------
void CalculateDamage::update()
{
    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();
    mat & data = m_particles->data();

    // Updating single particle states
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        const vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
        const int jnum = PDconnections.size();
        const double maxConnections = jnum;

        if(maxConnections <= 0)
        {
            data(i, m_indexDamage) = 0;
            return;
        }

        double tot = 0;
        for (const auto& con:PDconnections)
        {
            if(con.second[m_indexConnected] > 0.5)
                tot += 1;
        }

        data(i, m_indexDamage) = 1. - tot/maxConnections;
        //    data(i, m_indexDamage) = tot;
    }
}
//------------------------------------------------------------------------------
}
