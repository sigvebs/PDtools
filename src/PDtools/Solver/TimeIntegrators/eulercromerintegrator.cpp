#include "eulercromerintegrator.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
EulerCromerIntegrator::EulerCromerIntegrator() {}
//------------------------------------------------------------------------------
EulerCromerIntegrator::~EulerCromerIntegrator() {}
//------------------------------------------------------------------------------
void EulerCromerIntegrator::integrateStepOne() {
  mat &r = m_particles->r();
  mat &v = m_particles->v();
  const mat &F = m_particles->F();
  const mat &data = m_particles->data();
  const arma::imat &isStatic = m_particles->isStatic();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < m_particles->nParticles(); i++) {
    const int col_i = i;
    if (isStatic(col_i)) {
      continue;
    }
    const double rho = data(col_i, m_indexRho);
    const double dtRho = m_dt / rho;

    for (int d = 0; d < m_dim; d++) {
      const double Fd = F(col_i, d);
      v(col_i, d) += Fd * dtRho;
      r(col_i, d) += v(col_i, d) * m_dt;
    }
  }
}
//------------------------------------------------------------------------------
void EulerCromerIntegrator::integrateStepTwo() {}
//------------------------------------------------------------------------------
}
