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
    mat & F = m_particles->F();
    const mat & data = m_particles->data();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double rho = data(col_i, m_colRho);
        double dtRho = 0.5*m_dt/rho;

        for(int d=0; d<m_dim; d++)
        {
            v(d, col_i) += F(d, col_i)*dtRho;
            r(d, col_i) += v(d, col_i)*m_dt;
            F(d, col_i) = 0;
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

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double rho = data(col_i, m_colRho);
        double dtRho = 0.5*m_dt/rho;

        for(int d=0; d<m_dim; d++)
        {
            v(d, col_i) += F(d, col_i)*dtRho;
        }
    }
}
//------------------------------------------------------------------------------
}

