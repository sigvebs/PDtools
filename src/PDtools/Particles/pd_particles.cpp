#include "pd_particles.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
PD_Particles::PD_Particles()
{

}
//------------------------------------------------------------------------------
PD_Particles::~PD_Particles()
{

}
//------------------------------------------------------------------------------
void PD_Particles::initializeMatrices()
{
    Particles::initializeMatrices();

    m_r0 = mat(m_nParticles*PARTICLE_BUFFER, DIM);
    m_v  = mat(m_nParticles*PARTICLE_BUFFER, DIM);
    m_F  = mat(m_nParticles*PARTICLE_BUFFER, DIM);
}
//------------------------------------------------------------------------------
void PD_Particles::initializeADR()
{
    m_stableMass = vec(m_nParticles*PARTICLE_BUFFER);
//    m_Fold  = mat(DIM, m_nParticles*PARTICLE_BUFFER);
    m_Fold  = mat(m_nParticles*PARTICLE_BUFFER, DIM);
}
//------------------------------------------------------------------------------
void PD_Particles::initializeBodyForces()
{
    m_u = mat(m_nParticles, DIM); // TMP
    m_b = mat(m_nParticles, DIM); // TMP
//    m_b = mat(m_nParticles*PARTICLE_BUFFER, DIM);
}
//------------------------------------------------------------------------------
void PD_Particles::dimensionalScaling(const double E0, const double L0, const double v0,
                                      const double t0, const double rho0)
{
    (void) E0;
    (void) t0;

    vector<pair<int, double>> parameterScaling;

    if(hasParameter("rho"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("rho"), rho0));
    }
    if(hasParameter("mass"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("mass"), rho0));
    }
    if(hasParameter("volume"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("volume"), L0*L0*L0));
    }

    for(unsigned int i=0; i<m_nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            m_r(i, d) /= L0;
            m_r0(i, d) /= L0;
            m_v(i, d) /= v0;
        }

        for(auto pos_scaling:parameterScaling)
        {
            m_data(i, pos_scaling.first) /= pos_scaling.second;
        }
    }
}
//------------------------------------------------------------------------------
}
