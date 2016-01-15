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
void ComputeKineticEnergy::update(const int id_i, const int i)
{
    double rho = (*m_data)(i, m_indexRho);
    double volume = (*m_data)(i, m_indexVolume);
    double mass = rho*volume;

    double v_squared = 0;
    for(int d=0; d<m_dim; d++)
    {
        v_squared += (*m_v)(i, d)*(*m_v)(i, d);
    }

    (*m_data)(i, m_indexKE) = 0.5*mass*v_squared;
}
//------------------------------------------------------------------------------
}
