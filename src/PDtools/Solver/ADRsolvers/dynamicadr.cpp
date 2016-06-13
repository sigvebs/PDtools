#include "dynamicadr.h"

//#include "PDtools/Particles/pd_particles.h"
//#include "PDtools/Force/force.h"
//#include "PDtools/Modfiers/modifier.h"
//#include "PDtools/SavePdData/savepddata.h"
using namespace arma;

namespace PDtools
//------------------------------------------------------------------------------
{
dynamicADR::dynamicADR()
{

}
//------------------------------------------------------------------------------
void dynamicADR::solve()
{
    initialize();
    checkInitialization();
    calculateForces(0);
    updateProperties(0);
    save(0);

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }

}
//------------------------------------------------------------------------------
void dynamicADR::stepForward(int timeStep)
{
    integrateStepOne();
    updateGridAndCommunication();

    modifiersStepOne();
    updateGridAndCommunication();

    updateProperties(timeStep + 1);
    updateGridAndCommunication();

    save(timeStep + 1);
    //----------------------------------------------------------------------
    zeroForces();
    calculateForces(timeStep+1);
    //----------------------------------------------------------------------

    modifiersStepTwo();
    integrateStepTwo();

    m_t += m_dt;
}
//------------------------------------------------------------------------------
}


