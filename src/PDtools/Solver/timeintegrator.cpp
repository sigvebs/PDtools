#include "timeintegrator.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/SavePdData/savepddata.h"
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
void TimeIntegrator::solve()
{
    initialize();
    checkInitialization();

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void TimeIntegrator::stepForward(int i)
{
    integrateStepOne();
    updateGridAndCommunication();

    modifiersStepOne();

    calculateForces();

    modifiersStepTwo();
    integrateStepTwo();

    save(i);
    m_t += m_dt;

    //    modifiersStepOne();
//    integrateStepOne();

//    updateGridAndCommunication();

//    calculateForces();

//    integrateStepTwo();
//    modifiersStepTwo();

//    save(i);
//    m_t += m_dt;
}
//------------------------------------------------------------------------------
void TimeIntegrator::initialize()
{
    if(m_particles->hasParameter("rho"))
    {
        m_colRho = m_particles->getParamId("rho");
    }
    else if(m_particles->hasParameter("mass"))
    {
        m_colRho = m_particles->getParamId("mass");
    }
    else
    {
        cerr << "ERROR: Particles data does not contain either"
             << " 'rho' or 'mass'. This is needed for time integration." << endl;
        throw ParticlesMissingRhoOrMass;
    }

    Solver::initialize();
}
//------------------------------------------------------------------------------
void TimeIntegrator::checkInitialization()
{
    if(m_dt == 0)
    {
        cerr << "Error: 'dt' is not set in the timeintegrator" << endl;
        throw TimeStepNotSet;
    }

    // Check for forces and order of modifiers
    Solver::checkInitialization();
}
//------------------------------------------------------------------------------
}
