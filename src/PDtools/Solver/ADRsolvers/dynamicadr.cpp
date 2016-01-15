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

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void dynamicADR::stepForward(int timeStep)
{
    save(timeStep);
    modifiersStepOne();
    integrateStepOne();
    zeroForces();

    updateGridAndCommunication();

    calculateForces(timeStep);

    modifiersStepTwo();
    integrateStepTwo();

    m_t += m_dt;
    cout << timeStep << " dt = " << m_t << endl;
}
//------------------------------------------------------------------------------
}


