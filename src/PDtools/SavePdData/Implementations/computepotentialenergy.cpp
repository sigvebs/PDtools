#include "computepotentialenergy.h"

#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ComputePotentialEnergy::ComputePotentialEnergy(PD_Particles &particles, vector<Force *> &forces):
    ComputeProperty(particles),
    m_forces(forces),
    m_data(m_particles.data())
{
    m_indexPotential = m_particles.registerParameter("potential_energy");
}
//------------------------------------------------------------------------------
ComputePotentialEnergy::~ComputePotentialEnergy()
{

}
//------------------------------------------------------------------------------
void ComputePotentialEnergy::update(const int id_i, const int i)
{
    m_data(id_i, m_indexPotential) = 0;

    for(Force *force: m_forces)
    {
        force->calculatePotentialEnergy(id_i, i, m_indexPotential);
    }
}
//------------------------------------------------------------------------------
}

