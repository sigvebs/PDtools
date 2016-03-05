#include "solver.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/SavePdData/savepddata.h"
#include "PDtools/Modfiers/modifier.h"
#include "PDtools/PdFunctions/pdfunctions.h"
#include "PDtools/CalculateProperties/calculateproperty.h"
#include <iostream>

#include <armadillo>
using namespace arma;

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
void Solver::setRankAndCores(int rank, int cores)
{
    m_myRank = rank;
    m_nCores = cores;
}
//------------------------------------------------------------------------------
void Solver::setCalculateProperties(vector<CalculateProperty*> &calcProp)
{
    m_properties = calcProp;
}
//------------------------------------------------------------------------------
void Solver::setSaveParticles(SavePdData *saveParticles)
{
    m_saveParticles = saveParticles;
}
//------------------------------------------------------------------------------
Solver::Solver()
{
    
}
//------------------------------------------------------------------------------
Solver::~Solver()
{     
    for(Force *force:m_oneBodyForces)
    {
      delete force;
    } 
    m_oneBodyForces.clear();
}
//------------------------------------------------------------------------------
void Solver::updateGridAndCommunication()
{
    m_mainGrid->clearParticles();
    updateGrid(*m_mainGrid, *m_particles);

#if USE_MPI
    exchangeGhostParticles(*m_mainGrid, *m_particles);

    // Updating lists in modifiers
    int counter = 0;
    for(Modifier *modifier:m_modifiers)
    {
        updateModifierLists(*modifier, *m_particles, counter);
        counter++;
    }
#endif
}
//------------------------------------------------------------------------------
void Solver::updateGhosts()
{
    // TODO: needs optimization
    m_mainGrid->clearGhostParticles();
    exchangeGhostParticles(*m_mainGrid, *m_particles);
}
//------------------------------------------------------------------------------
void Solver::save(int timesStep)
{
    if(timesStep % m_saveInterval == 0)
    {
        m_saveParticles->evaluate(m_t, timesStep);
        m_saveParticles->saveData(m_t, timesStep);

        const double progress = (double)timesStep/(double)m_steps;
        if(m_myRank == 0)
            printProgress(progress);
    }
}
//------------------------------------------------------------------------------
void Solver::initialize()
{
}
//------------------------------------------------------------------------------
void Solver::modifiersStepOne()
{
    for(Modifier *modifier:m_modifiers)
    {
        modifier->evaluateStepOne();
    }

    for(Modifier *modifier:m_spModifiers)
    {
        modifier->evaluateStepOne();
    }
    const ivec &colToId = m_particles->colToId();
    const int nParticles =  m_particles->nParticles();
    if(!m_spModifiers.empty())
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<nParticles; i++)
        {
            const int id = colToId(i);
            for(Modifier *modifier:m_spModifiers)
            {
                modifier->evaluateStepOne(id, i);
            }
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<nParticles; i++)
        {
            const int id = colToId(i);
            for(Modifier *modifier:m_spModifiers)
            {
                modifier->updateStepOne(id, i);
            }
        }
    }

    for(Modifier *modifier:m_spModifiers)
    {
        modifier->evaluateStepOnePost();
    }
}
//------------------------------------------------------------------------------
void Solver::modifiersStepTwo()
{
    for(Modifier *modifier:m_modifiers)
    {
        modifier->evaluateStepTwo();
    }

    const ivec &colToId = m_particles->colToId();
    const int nParticles =  m_particles->nParticles();

    if(!m_spModifiers.empty())
    {
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<nParticles; i++)
        {
            const int id = colToId(i);
            for(Modifier *modifier:m_spModifiers)
            {
                modifier->evaluateStepTwo(id, i);
            }
        }
    }
}
//------------------------------------------------------------------------------
void Solver::zeroForces()
{
    const int nParticles = m_particles->nParticles() + m_particles->nGhostParticles();
    mat & F = m_particles->F();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            F(i, d) = 0;
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
void Solver::setDim(double _dim)
{
    m_dim = _dim;
}
//------------------------------------------------------------------------------
void Solver::addForce(Force *force)
{
    m_oneBodyForces.push_back(force);
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
void Solver::addQsModifiers(Modifier *modifier)
{
    m_qsModifiers.push_back(modifier);
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
void Solver::calculateForces(int timeStep)
{
    const ivec &colToId = m_particles->colToId();
    const imat & isStatic = m_particles->isStatic();
    const int nParticles = m_particles->nParticles();

    // Updating overall state
    for(Force *oneBodyForce:m_oneBodyForces)
    {
        oneBodyForce->updateState();
    }
    // Updating single particle states
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id = colToId(i);
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            oneBodyForce->updateState(id, i);
        }
    }

    // Calculating one-body forces
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id = colToId(i);

        if(isStatic(i))
            continue;

        for(Force *oneBodyForce:m_oneBodyForces)
        {
            oneBodyForce->calculateForces(id, i);
        }
    }
}
//------------------------------------------------------------------------------
void Solver::updateProperties(const int timeStep)
{
    for(CalculateProperty* property:m_properties)
    {
        const int updateFrequency = property->updateFrquency();

        if(timeStep % updateFrequency == 0)
        {
            property->clean();
            property->update();
        }
    }
}
//------------------------------------------------------------------------------
void Solver::printProgress(const double progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos)
            std::cout << "=";
        else if (j == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
//------------------------------------------------------------------------------
}
