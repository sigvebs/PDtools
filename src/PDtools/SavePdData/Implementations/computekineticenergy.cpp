#include "computekineticenergy.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ComputeKineticEnergy::ComputeKineticEnergy(PD_Particles &particles):
    ComputeProperty(particles)
{
    m_indexKE = m_particles.registerParameter("kinetic_energy");
    m_indexVolume =  m_particles.getParamId("volume");
    m_indexRho =  m_particles.getParamId("rho");
    m_v = &m_particles.v();
    m_data = &m_particles.data();
}
//------------------------------------------------------------------------------
ComputeKineticEnergy::~ComputeKineticEnergy()
{

}

//------------------------------------------------------------------------------
void ComputeKineticEnergy::update(const pair<int, int> &pIdcol)
{
    int col = pIdcol.second;
    double rho = (*m_data)(col, m_indexRho);
    double volume = (*m_data)(col, m_indexVolume);
    double mass = rho*volume;

    double v_squared = 0;
    for(int d=0; d<m_dim; d++)
    {
        v_squared += (*m_v)(d, col)*(*m_v)(d, col);
    }

    (*m_data)(col, m_indexKE) = 0.5*mass*v_squared;
}
//------------------------------------------------------------------------------
}
