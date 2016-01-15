#include "adr.h"

#if USE_MPI
#include <mpi.h>
#endif

#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/Modfiers/modifier.h"
#include "PDtools/SavePdData/savepddata.h"
#include "PDtools/Grid/grid.h"
#include "PDtools/PdFunctions/pdfunctions.h"

using namespace arma;

namespace PDtools
//------------------------------------------------------------------------------
{
ADR::ADR()
{
    m_dt = 1.0;
    m_c = 0;
}
//------------------------------------------------------------------------------
ADR::~ADR()
{

}
//------------------------------------------------------------------------------
void ADR::solve()
{
    initialize();
    checkInitialization();
    save(0);
    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void ADR::stepForward(int i)
{
//    cerr << "void ADR::stepForward(int i): must be update with calc properties" << endl;
//    exit(1);
    modifiersStepOne();

    m_globalError = 2*m_errorThreshold;
    iterate(m_maxSteps);
    updateProperties(i);

    if(m_myRank == 0)
        cout << "i = " << i << " " << m_globalError << endl;
    modifiersStepTwo();
    save(i+1);
    m_t += m_dt;
}
//------------------------------------------------------------------------------
void ADR::iterate(int maxNumberOfSteps)
{
    int counter = 0;
    int minSteps = 5;
    int nFractureRelaxationSteps = m_maxStepsFracture;
    const ivec &colToId = m_particles->colToId();

    do
    {
        integrateStepOne();
        zeroForces();
        updateGridAndCommunication();
        calculateForces(0);
        staticModifiers();
        integrateStepTwo();

        counter++;
        if(counter > maxNumberOfSteps)
            break;
    }
    while( (m_globalError > m_errorThreshold) || (counter < minSteps));
    if(m_myRank == 0)
        cout << "counter:" << counter << endl;
    int nParticles;
#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //--------------------------------------------------
    if(!m_qsModifiers.empty())
    {
        nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<nParticles; i++)
        {
            const int id = colToId(i);
            for(Modifier *modifier:m_qsModifiers)
            {
                modifier->evaluateStepOne(id, i);
            }
        }
        updateGridAndCommunication();
        updateProperties(0);
        updateGridAndCommunication();
        nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        for(int i=0; i<nParticles; i++)
        {
            const int id = colToId(i);
            for(Modifier *modifier:m_qsModifiers)
            {
                modifier->evaluateStepTwo(id, i);
            }

        }
        for(Modifier *modifier:m_qsModifiers)
        {
            modifier->evaluateStepTwo();
        }

        int continueState = false;

        for(Modifier *modifier:m_qsModifiers)
        {
            if(modifier->state())
                continueState = 1;
        }

#if USE_MPI
    MPI_Allreduce(&continueState, &continueState, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
        if(continueState > 0)
        {
            m_globalError = 2*m_errorThreshold;
            iterate(nFractureRelaxationSteps);
        }
        else
        {
            if(m_myRank == 0)
               cout << "counter = " << counter << " error = " << m_globalError << endl;
        }
    }
//    int counter = 0;
//    int minSteps = 5;
//    int nFractureRelaxationSteps = m_maxStepsFracture;
//    const ivec &colToId = m_particles->colToId();

//    do
//    {
//        integrateStepOne();
//        zeroForcesAndStress();

//        updateGridAndCommunication();
//        calculateForces(0);
//        staticModifiers();

//        integrateStepTwo();

//        counter++;
//        if(counter > maxNumberOfSteps)
//            break;
//    }
//    while( (m_globalError > m_errorThreshold) || (counter < minSteps));

//    const int nParticles = m_particles->nParticles();
//#if USE_MPI
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif
//    //--------------------------------------------------
//    if(!m_qsModifiers.empty())
//    {
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
//        for(int i=0; i<nParticles; i++)
//        {
//            const int id = colToId(i);
//            for(Modifier *modifier:m_qsModifiers)
//            {
//                modifier->evaluateStepOne(id, i);
//            }
//        }
//        updateGridAndCommunication();
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
//        for(int i=0; i<nParticles; i++)
//        {
//            const int id = colToId(i);
//            for(Modifier *modifier:m_qsModifiers)
//            {
//                modifier->evaluateStepTwo(id, i);
//            }

//        }
//        for(Modifier *modifier:m_qsModifiers)
//        {
//            modifier->evaluateStepTwo();
//        }

//        int continueState = false;

//        for(Modifier *modifier:m_qsModifiers)
//        {
//            if(modifier->state())
//                continueState = 1;
//        }

//#if USE_MPI
//    MPI_Allreduce(&continueState, &continueState, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//#endif
//        if(continueState > 0)
//        {
//            m_globalError = 2*m_errorThreshold;
//            iterate(nFractureRelaxationSteps);
//        }
//        else
//        {
//            if(m_myRank == 0)
//               cout << "counter = " << counter << " error = " << m_globalError << endl;
//        }
//    }
}
//------------------------------------------------------------------------------
void ADR::checkInitialization()
{
}
//------------------------------------------------------------------------------
void ADR::initialize()
{
    //-----------------------------
    mat & F = m_particles->F();
    mat & v = m_particles->v();
    mat & Fold = m_particles->Fold();
    int nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        for(int d=0; d<m_dim; d++){
            F(i, d) = 0;
            Fold(i, d) = 0;
            v(i, d) = 0;
        }
    }

    updateGridAndCommunication();
    //-----------------------------
    calculateForces(0);
    calculateStableMass();
    Solver::initialize();
}
//------------------------------------------------------------------------------
void ADR::calculateStableMass()
{
    const ivec &colToId = m_particles->colToId();
    arma::vec &stableMass = m_particles->stableMass();
    const int nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id = colToId(i);
        double stableMassMax = 0;
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            const double sm = oneBodyForce->calculateStableMass(id, i, m_dt);
            stableMassMax = max(stableMassMax, sm);
        }
        stableMass(i) = stableMassMax;
    }
}
//------------------------------------------------------------------------------
void ADR::integrateStepOne()
{
    using namespace arma;

    const double alpha = (2.*m_dt)/(2. + m_c*m_dt);
    const double beta = (2. - m_c*m_dt)/(2. + m_c*m_dt);

    const arma::imat & isStatic = m_particles->isStatic();
    const vec & stableMass = m_particles->stableMass();
    mat & r = m_particles->r();
    mat & v = m_particles->v();
    mat & F = m_particles->F();
    mat & Fold = m_particles->Fold();

    double Fn = 0;
    double deltaFn = 0;
    const int nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+: deltaFn, Fn)
#endif
    for(int i=0; i<nParticles; i++)
    {
        if(isStatic(i))
        {
            continue;
        }
        for(int d=0; d<m_dim; d++){
            v(i, d) = beta*v(i, d) + alpha*F(i, d)/stableMass(i);
            r(i, d) += v(i, d)*m_dt;

            deltaFn += pow(Fold(i, d) - F(i, d), 2);
            Fn += pow(F(i, d), 2);

            Fold(i, d) = F(i, d);
        }
    }
#if USE_MPI
    MPI_Allreduce(&deltaFn, &deltaFn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Fn, &Fn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    m_globalError = sqrt(deltaFn/Fn);

    if(std::isnan(m_globalError))
        m_globalError = 2*m_errorThreshold;
}
//------------------------------------------------------------------------------
void ADR::integrateStepTwo()
{
    const mat & v = m_particles->v();
    const mat & F = m_particles->F();
    const mat & Fold = m_particles->Fold();
    const vec & stableMass = m_particles->stableMass();
    const arma::imat & isStatic = m_particles->isStatic();

    // Calculating the damping coefficient
    double numerator = 0;
    double denominator = 0;
    int const nParticles = m_particles->nParticles();

#ifdef USE_OPENMP
# pragma omp parallel for reduction(+: numerator, denominator)
#endif
    for(int i=0; i<nParticles; i++)
    {
        if(isStatic(i))
            continue;
        for(int d=0; d<m_dim; d++)
        {
            numerator += -v(i, d)*(F(i, d) - Fold(i, d))
                    /(stableMass(i)*m_dt);
            denominator += v(i, d)*v(i, d);
        }
    }

#if USE_MPI
    MPI_Allreduce(&numerator, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&denominator, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    m_c = 0;
    if(denominator > 0)
    {
        if(numerator/denominator > 0)
        {
            m_c = 2*sqrt(numerator/denominator);
        }
    }

//    if(m_c >= 2.0)
//    {
//        m_c = 1.99;
//    }
}
//------------------------------------------------------------------------------
void ADR::staticModifiers()
{
    for(Modifier *modifier:m_modifiers)
    {
        modifier->staticEvaluation();
    }
}
//------------------------------------------------------------------------------
void ADR::updateGridAndCommunication()
{
    m_mainGrid->clearParticles();
    updateGrid(*m_mainGrid, *m_particles, true);

#if USE_MPI
    m_mainGrid->clearGhostParticles();
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
}

