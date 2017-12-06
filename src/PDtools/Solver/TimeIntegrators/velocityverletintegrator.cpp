#include "velocityverletintegrator.h"

#include <armadillo>
using namespace arma;

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
VelocityVerletIntegrator::VelocityVerletIntegrator()
{

}
//------------------------------------------------------------------------------
VelocityVerletIntegrator::~VelocityVerletIntegrator()
{

}
//------------------------------------------------------------------------------
void VelocityVerletIntegrator::integrateStepOne()
{
    // Calulating the half step velocity for all particeles
    // v(t + 0.5dt) = v(t) + 0.5 (f(t)V + b(t))/rho dt
    // x(t + dt)    = x(t) + v(t + 0.5dt) dt
    mat & r = m_particles->r();
    mat & v = m_particles->v();
    const mat & F = m_particles->F();
    const mat & data = m_particles->data();
    const double dtRhoHalf = 0.5*m_dt;
    const arma::imat & isStatic = m_particles->isStatic();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++) {
        if(isStatic(i)) {
            continue;
        }
        const double rho = data(i, m_indexRho);
        const double dtRho = dtRhoHalf/rho;

        for(int d=0; d<m_dim; d++) {
            double Fd = F(i, d);
            Fd *= dtRho;
            v(i, d) += Fd;
            r(i, d) += v(i, d)*m_dt;
        }
    }
}
//------------------------------------------------------------------------------
void VelocityVerletIntegrator::integrateStepTwo()
{
    // Updating to the full-step velocity
    mat & v = m_particles->v();
    const mat & F = m_particles->F();
    const mat & data = m_particles->data();
    const double dtRhoHalf = 0.5*m_dt;
    const arma::imat & isStatic = m_particles->isStatic();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++) {
        if(isStatic(i))
            continue;

        const double rho = data(i, m_indexRho);
        const double dtRho = dtRhoHalf/rho;

        for(int d=0; d<m_dim; d++) {
            v(i, d) += F(i, d)*dtRho;
        }
    }
}
//------------------------------------------------------------------------------
}

