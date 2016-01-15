#include "modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
std::vector<pair<string, int> > Modifier::neededProperties() const
{
    return m_neededProperties;
}
//------------------------------------------------------------------------------
int Modifier::dim() const
{
    return m_dim;
}
//------------------------------------------------------------------------------
void Modifier::setDim(int dim)
{
    m_dim = dim;
}
//------------------------------------------------------------------------------
Modifier::Modifier()
{

}
//------------------------------------------------------------------------------
Modifier::~Modifier()
{

}
//------------------------------------------------------------------------------
void Modifier::registerParticleParameters()
{

}
//------------------------------------------------------------------------------
void Modifier::evaluateStepOne(const int id, const int i)
{
    (void) id;
    (void) i;
    return;
}
//------------------------------------------------------------------------------
void Modifier::updateStepOne(const int id, const int i)
{
    (void) id;
    (void) i;
}
//------------------------------------------------------------------------------
void Modifier::evaluateStepTwo(const int id, const int i)
{
    (void) id;
    (void) i;
    return;
}
//------------------------------------------------------------------------------
void Modifier::evaluateStepOne()
{

}
//------------------------------------------------------------------------------
void Modifier::evaluateStepTwo()
{

}
//------------------------------------------------------------------------------
void Modifier::staticEvaluation()
{
    return;
}
//------------------------------------------------------------------------------
void Modifier::initialize()
{
    return;
}

void Modifier::setParticles(PD_Particles &particles)
{
    m_particles = &particles;
}
//------------------------------------------------------------------------------
bool Modifier::state()
{
    return m_state;
}
//------------------------------------------------------------------------------
const std::vector<string> &Modifier::initalGhostDependencies()
{
    return m_initialGhostParameters;
}
//------------------------------------------------------------------------------
const std::vector<string> &Modifier::ghostDependencies()
{
    return m_ghostParameters;
}
//------------------------------------------------------------------------------
void Modifier::addToList(int id)
{
    m_localParticleIds.push_back(id);
}
//------------------------------------------------------------------------------
bool Modifier::removeFromList(const int id)
{
    bool found = false;
    int counter = -1;

    for(const int &l_id:m_localParticleIds)
    {
        counter++;
        if(l_id == id)
        {
            found = true;
            break;
        }
    }
    if(found)
    {
        m_localParticleIds.erase(m_localParticleIds.begin() + counter);
    }

    return found;
}
//------------------------------------------------------------------------------
}
