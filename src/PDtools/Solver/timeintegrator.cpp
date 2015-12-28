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
    calculateForces();

    cout << "Starting time integration " << endl;

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void TimeIntegrator::stepForward(int i)
{
    save(i);
    modifiersStepOne();
    integrateStepOne();
    zeroForcesAndStress();

    updateGridAndCommunication();

    calculateForces();

    modifiersStepTwo();
    integrateStepTwo();

    m_t += m_dt;
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
