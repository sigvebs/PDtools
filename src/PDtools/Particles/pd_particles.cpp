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

    m_r0 = arma::mat(DIM, m_nParticles*PARTICLE_BUFFER);
    m_v  = arma::mat(DIM, m_nParticles*PARTICLE_BUFFER);
    m_F  = arma::mat(DIM, m_nParticles*PARTICLE_BUFFER);
}
//------------------------------------------------------------------------------
void PD_Particles::initializeADR()
{
    m_stableMass = arma::vec(m_nParticles*PARTICLE_BUFFER);
    m_Fold  = arma::mat(DIM, m_nParticles*PARTICLE_BUFFER);
}
//------------------------------------------------------------------------------
}
