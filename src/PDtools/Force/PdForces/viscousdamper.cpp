#include "viscousdamper.h"

namespace PDtools {
//------------------------------------------------------------------------------
ViscousDamper::ViscousDamper(PDtools::PD_Particles &particles,
                             const double dampeningCoeff)
    : Force(particles), m_c(dampeningCoeff), m_v(m_particles.v()),
      m_isStatic(m_particles.isStatic()) {}
//------------------------------------------------------------------------------
void ViscousDamper::calculateForces(const int id_i, const int i) {
  (void)id_i;

  if (m_isStatic(id_i))
    return;

  for (int d = 0; d < m_dim; d++) {
    m_v(i, d) *= m_c;
  }
}
//------------------------------------------------------------------------------
}
