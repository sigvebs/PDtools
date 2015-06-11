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
}
//------------------------------------------------------------------------------
ComputeDamage::~ComputeDamage()
{

}
//------------------------------------------------------------------------------
void ComputeDamage::update(const pair<int, int> &pIdcol)
{
    double maxConnections = (*m_data)(pIdcol.second, m_indexMaxPdConnections);
    (*m_data)(pIdcol.second, m_indexDamage) = 1. - m_particles.pdConnections(pIdcol.first).size()/maxConnections;
}
//------------------------------------------------------------------------------
void ComputeDamage::init(const pair<int, int> &pIdcol)
{
    (*m_data)(pIdcol.second, m_indexMaxPdConnections) = m_particles.pdConnections(pIdcol.first).size();
}
//------------------------------------------------------------------------------

}

