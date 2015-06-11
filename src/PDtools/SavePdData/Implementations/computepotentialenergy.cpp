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
void ComputePotentialEnergy::update(const pair<int, int> &pIdcol)
{
    m_data(pIdcol.second, m_indexPotential) = 0;

    for(Force *force: m_forces)
    {
        force->calculatePotentialEnergy(pIdcol, m_indexPotential);
    }
}
//------------------------------------------------------------------------------
}

