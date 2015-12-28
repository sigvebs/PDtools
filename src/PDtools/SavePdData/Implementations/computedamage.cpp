#include "computedamage.h"


namespace PDtools
{
//------------------------------------------------------------------------------
ComputeDamage::ComputeDamage(PD_Particles &particles):
    ComputeProperty(particles)
{
    m_indexDamage = m_particles.registerParameter("damage");
    m_indexMaxPdConnections = m_particles.registerParameter("maxPdConnections");
    m_data = &m_particles.data();
    m_indexConnected = m_particles.getPdParamId("connected");
}
//------------------------------------------------------------------------------
ComputeDamage::~ComputeDamage()
{

}
//------------------------------------------------------------------------------
void ComputeDamage::update(const pair<int, int> &pIdcol)
{
    const double maxConnections = (*m_data)(pIdcol.second, m_indexMaxPdConnections);
    if(maxConnections <= 0)
    {
        (*m_data)(pIdcol.second, m_indexDamage) = 0;
        return;
    }

    const int pId = pIdcol.first;
    const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    const int jnum = PDconnections.size();

    double tot = 0;
    for (int jj = 0; jj < jnum; jj++)
    {
        const auto &con = PDconnections[jj];
        if(con.second[m_indexConnected] <= 0.5)
            continue;
        tot += 1;
    }

    (*m_data)(pIdcol.second, m_indexDamage) = 1. - tot/maxConnections;
}
//------------------------------------------------------------------------------
void ComputeDamage::init(const pair<int, int> &pIdcol)
{
    (*m_data)(pIdcol.second, m_indexMaxPdConnections) = m_particles.pdConnections(pIdcol.first).size();
}
//------------------------------------------------------------------------------

}

