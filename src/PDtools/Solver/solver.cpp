#include "solver.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/SavePdData/savepddata.h"
#include "PDtools/Modfiers/modifier.h"
#include <iostream>

namespace PDtools
{
//------------------------------------------------------------------------------
void Solver::setADR_fracture(Modifier *ADR_fracture)
{
    m_ADR_fracture = ADR_fracture;
}
//------------------------------------------------------------------------------
void Solver::setErrorThreshold(double errorThreshold)
{
    m_errorThreshold = errorThreshold;
}
//------------------------------------------------------------------------------
Solver::Solver()
{

}
//------------------------------------------------------------------------------
Solver::~Solver()
{

}
//------------------------------------------------------------------------------
std::vector<std::string> Solver::saveParameters() const
{
    return m_saveParameters;
}
//------------------------------------------------------------------------------
void Solver::setSaveParameters(const std::vector<std::string> &saveParameters)
{
    m_saveParameters = saveParameters;
}
//------------------------------------------------------------------------------
void Solver::updateGridAndCommunication()
{
    m_mainGrid->placeParticlesInGrid(*m_particles);
    // MPI COMMUNICATIONS HERE
}
//------------------------------------------------------------------------------
void Solver::save(int i)
{
    if(i%m_saveInterval == 0)
    {
        saveParticles->evaluate(m_t, i);
    }
}
//------------------------------------------------------------------------------
void Solver::initialize()
{
    saveParticles = new SavePdData(m_saveParameters);
    saveParticles->setSavePath(m_savePath);
    saveParticles->setParticles(m_particles);
    saveParticles->setForces(m_oneBodyForces);
    saveParticles->initialize();
}
//------------------------------------------------------------------------------
void Solver::modifiersStepOne()
{
    for(Modifier *modifier:m_modifiers)
    {
        modifier->evaluateStepOne();
    }

    if(!m_spModifiers.empty())
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<m_particles->nParticles(); i++)
        {
            pair<int, int> id(i, i);
            for(Modifier *modifier:m_spModifiers)
            {
                modifier->evaluateStepOne(id);
            }

        }
    }
}
//------------------------------------------------------------------------------
void Solver::modifiersStepTwo()
{
    for(Modifier *modifier:m_modifiers)
    {
        modifier->evaluateStepTwo();
    }

    if(!m_spModifiers.empty())
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<m_particles->nParticles(); i++)
        {
            pair<int, int> id(i, i);
            for(Modifier *modifier:m_spModifiers)
            {
                modifier->evaluateStepTwo(id);
            }
        }
    }
}
//------------------------------------------------------------------------------
void Solver::setParticles(PD_Particles &_particles)
{
    m_particles = &_particles;
}
//------------------------------------------------------------------------------
void Solver::setDomain(Domain &_domain)
{
    m_domain = &_domain;
}
//------------------------------------------------------------------------------
void Solver::setMainGrid(Grid &_grid)
{
    m_mainGrid = &_grid;
}
//------------------------------------------------------------------------------
void Solver::setSteps(int steps)
{
    m_steps = steps;
}
//------------------------------------------------------------------------------
void Solver::setDt(double _dt)
{
    m_dt = _dt;
}
//------------------------------------------------------------------------------
void Solver::setT(double _t)
{
    m_t = _t;
}
//------------------------------------------------------------------------------
void Solver::addForce(Force *force)
{
    m_oneBodyForces.push_back(force);
}
//------------------------------------------------------------------------------
void Solver::setSavePath(const string &savePath)
{
    m_savePath = savePath;
}
//------------------------------------------------------------------------------
void Solver::setSaveInterval(double saveInterval)
{
    m_saveInterval = saveInterval;
}
//------------------------------------------------------------------------------
void Solver::addSpModifier(Modifier * modifier)
{
    m_spModifiers.push_back(modifier);
}
//------------------------------------------------------------------------------
void Solver::addModifier(Modifier *modifier)
{
    m_modifiers.push_back(modifier);
}
//------------------------------------------------------------------------------
void Solver::checkInitialization()
{
    if(m_steps == 0)
    {
        cerr << "Error: 'steps' is not set in the solver" << endl;
        throw NumberOfStepNotSet;
    }

    if(m_particles == nullptr)
    {
        cerr << "Error: 'particles' is not set in the solver" << endl;
        throw ParticlesNotSet;
    }

    if(m_mainGrid == nullptr)
    {
        cerr << "Error: 'mainGrid' is not set in the solver" << endl;
        throw ParticlesNotSet;
    }
}
//------------------------------------------------------------------------------
void Solver::calculateForces()
{
    int nParticles = m_particles->nParticles();

    for(Force *oneBodyForce:m_oneBodyForces)
    {
        oneBodyForce->updateState();
    }

    // Calculate one-body forces
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        pair<int, int> id(i, i);
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            oneBodyForce->calculateForces(id);
        }
    }
}
//------------------------------------------------------------------------------
}
