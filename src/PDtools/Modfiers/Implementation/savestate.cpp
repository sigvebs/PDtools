#include "savestate.h"

#include "PDtools/Particles/pd_particles.h"
namespace PDtools
{
//------------------------------------------------------------------------------
SaveState::SaveState(int frequency):
    m_frequency(frequency)
{
    m_counter = 0;
}
//------------------------------------------------------------------------------
void SaveState::initialize()
{

}
//------------------------------------------------------------------------------
void SaveState::evaluateStepOne()
{
    if(m_counter % m_frequency == 0)
    {


    }

    m_counter++;
}
//------------------------------------------------------------------------------
}
