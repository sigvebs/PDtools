#include "viscousdamper.h"

#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ViscousDamper::ViscousDamper(PDtools::PD_Particles &particles,
                                      const double dampeningCoeff):
    Force(particles),
    m_c(dampeningCoeff),
    m_v(m_particles.v()),
    m_isStatic(m_particles.isStatic())
{
}
//------------------------------------------------------------------------------
ViscousDamper::~ViscousDamper()
{

}
//------------------------------------------------------------------------------
void ViscousDamper::calculateForces(const int id, const int i)
{
    (void) id;

    if(m_isStatic(id))
        return;

    for(int d=0; d<m_dim; d++)
    {
        m_v(i, d) *= m_c;
    }
}
//------------------------------------------------------------------------------
}
