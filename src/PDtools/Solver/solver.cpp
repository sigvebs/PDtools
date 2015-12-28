#include "solver.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/SavePdData/savepddata.h"
#include "PDtools/Modfiers/modifier.h"
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
        saveParticles->saveData(m_t, i);

        const double progress = (double)i/(double)m_steps;
        printProgress(progress);
    }
}
//------------------------------------------------------------------------------
void Solver::initialize()
{
    saveParticles = new SavePdData(m_saveParameters);
    saveParticles->setScaling(m_E0, m_L0, m_v0, m_t0, m_rho0);
    saveParticles->setSavePath(m_savePath);
    saveParticles->setParticles(m_particles);
    saveParticles->setForces(m_oneBodyForces);
    saveParticles->initialize();

    if(m_particles->hasParameter("s_xx"))
        m_indexStress[0] = m_particles->getParamId("s_xx");
    else
        m_indexStress[0] = m_particles->registerParameter("s_xx");

    if(m_particles->hasParameter("s_yy"))
        m_indexStress[1] = m_particles->getParamId("s_yy");
    else
        m_indexStress[1] = m_particles->registerParameter("s_yy");

    if(m_particles->hasParameter("s_zz"))
        m_indexStress[2] = m_particles->getParamId("s_zz");
    else
        m_indexStress[2] = m_particles->registerParameter("s_zz");

    if(m_particles->hasParameter("s_xy"))
        m_indexStress[3] = m_particles->getParamId("s_xy");
    else
        m_indexStress[3] = m_particles->registerParameter("s_xy");

    if(m_particles->hasParameter("s_xz"))
        m_indexStress[4] = m_particles->getParamId("s_xz");
    else
        m_indexStress[4] = m_particles->registerParameter("s_xz");

    if(m_particles->hasParameter("s_yz"))
        m_indexStress[5] = m_particles->getParamId("s_yz");
    else
        m_indexStress[5] = m_particles->registerParameter("s_yz");
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
        for(unsigned int i=0; i<m_particles->nParticles(); i++)
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
        for(unsigned int i=0; i<m_particles->nParticles(); i++)
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
void Solver::zeroForcesAndStress()
{
    mat & F = m_particles->F();
    mat & data = m_particles->data();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        pair<int, int> id(i, i);
        for(int d=0; d<m_dim; d++)
        {
            F(i, d) = 0;
        }
        for(int s=0; s<6; s++)
        {
            data(i, m_indexStress[s]) = 0;
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
void Solver::setSaveScaling(const double E0, const double L0, const double v0, const double t0, const double rho0)
{
    m_E0 = E0;
    m_L0 = L0;
    m_v0 = v0;
    m_t0 = t0;
    m_rho0 = rho0;
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
void Solver::calculateForces()
{
    const imat & isStatic = m_particles->isStatic();
    int nParticles = m_particles->nParticles();

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
        pair<int, int> id(i, i);
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            oneBodyForce->updateState(id);
        }
    }

    // Calculateing the one-body forces
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        pair<int, int> id(i, i);
        if(isStatic(i))
            continue;
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            oneBodyForce->calculateForces(id);
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
