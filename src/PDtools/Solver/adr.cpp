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
    m_du_u = 0;
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
void ADR::checkInitialization()
{
}
//------------------------------------------------------------------------------
void ADR::initialize()
{
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
void ADR::stepForward(int i)
{
    modifiersStepOne();
    applyBoundaryConditions();
    updateGridAndCommunication();

    m_globalError = 2*m_errorThreshold;
    iterate(m_maxSteps);
    updateProperties(i);

    updateGridAndCommunication();
    modifiersStepTwo();
    updateGridAndCommunication();
    save(i+1);
    m_t += m_dt;

    if(m_myRank == 0)
        cout << "i = " << i << " " << m_globalError << endl;
}
//------------------------------------------------------------------------------
void ADR::iterate(int maxNumberOfSteps)
{
    //    int nFractureRelaxationSteps = m_maxStepsFracture;
    //    int nFractureRelaxationSteps = 2*minSteps;
    const ivec &colToId = m_particles->colToId();
    int minSteps = 20;
    int nRelaxSteps = 20;
    int counter = 0;
    int steps = 0;
    int relaxSteps = minSteps;

    for(counter=0; counter<maxNumberOfSteps; counter++)
    {
        integrateStepOne();
        updateGridAndCommunication();
        zeroForces();
        calculateForces(0);
        integrateStepTwo();
        staticModifiers();
        updateGridAndCommunication();

        if(relaxSteps == 0)
        {
            int continueState = 0;
            int nParticles = m_particles->nParticles();

            // It now loops over all properties, not only the ones needed in
            // static computations.
            updateProperties(0);
            updateGhosts();

            // Forces modifiers
            for(Force *oneBodyForce:m_oneBodyForces)
            {
                if(!oneBodyForce->getHasStaticModifier())
                    continue;

                for(int i=0; i<nParticles; i++)
                {
                    const int id = colToId(i);
                    oneBodyForce->evaluateStatic(id, i);
                }
                if(oneBodyForce->getContinueState() > 0 )
                    continueState = 1;
            }


            for(Modifier *modifier:m_qsModifiers)
            {
                if(!modifier->hasStepOne())
                    continue;

                for(int i=0; i<nParticles; i++)
                {
                    const int id = colToId(i);
                    modifier->evaluateStepOne(id, i);
                }
            }

            updateGhosts();
            nParticles = m_particles->nParticles();

            for(Modifier *modifier:m_qsModifiers)
            {
                if(!modifier->hasStepTwo())
                    continue;

                for(int i=0; i<nParticles; i++)
                {
                    const int id = colToId(i);
                    modifier->evaluateStepTwo(id, i);
                }

            }

            for(Modifier *modifier:m_qsModifiers)
            {
                modifier->evaluateStepTwo();
            }

            for(Modifier *modifier:m_qsModifiers)
            {
                if(modifier->state())
                    continueState = 1;
            }

#if USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, &continueState, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
            if(continueState > 0)
            {
                m_globalError = 2*m_errorThreshold;
                counter -= m_maxStepsFracture;
                relaxSteps = nRelaxSteps;
                if(counter < 0)
                    counter = 0;
            }
        }
        relaxSteps--;
        if(relaxSteps < 0)
            relaxSteps = 0;
        steps++;
        if( (m_globalError < m_errorThreshold) && (counter > minSteps))
            break;
    }

    if(m_myRank == 0)
        cout << "counter:" << steps << " err:" << m_globalError <<  endl;
}
//------------------------------------------------------------------------------
void ADR::integrateStepOne()
{
//    ivec & colToId = m_particles->colToId();
    using namespace arma;

    const double alpha = (2.*m_dt)/(2. + m_c*m_dt);
    const double beta = (2. - m_c*m_dt)/(2. + m_c*m_dt);

    const arma::imat & isStatic = m_particles->isStatic();
    const vec & stableMass = m_particles->stableMass();
    mat & r = m_particles->r();
    mat & v = m_particles->v();
    mat & F = m_particles->F();
    mat & Fold = m_particles->Fold();
    mat & r0 = m_particles->r0();

    const int nParticles = m_particles->nParticles();

    double maxU = 0;
    double avgU = 0;
    int np = 0;

#ifdef USE_OPENMP
# pragma omp parallel for reduction(+: deltaFn, Fn)
#endif
    for(int i=0; i<nParticles; i++)
    {
        if(isStatic(i))
        {
            continue;
        }

        double l_u = 0;

        for(int d=0; d<m_dim; d++)
        {
            v(i, d) = beta*v(i, d) + alpha*F(i, d)/stableMass(i);
            r(i, d) += v(i, d)*m_dt;

            l_u += pow(r(i, d) - r0(i,d), 2);

            Fold(i, d) = F(i, d);
        }

        if(l_u > 1e-22)
        {
            const double sqrt_lu = sqrt(l_u);
            maxU = max(sqrt_lu, maxU);
            avgU += sqrt_lu;
            np++;
        }
    }
#if USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &maxU, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &avgU, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

//    if(m_myRank == 0)
//        cout << maxU << "\t" <<  avgU << "\t" <<  maxU/avgU << "\t" << np << "\t" << m_globalError << endl;
    const int nTot = m_particles->totParticles();

    if (m_counter % 10 == 0 && np/nTot > 0.1)
    {
        avgU /= np;
        const double de = fabs(maxU/avgU - m_du_u);
        m_globalError = de;
        m_du_u = maxU/avgU;

        if(std::isnan(m_globalError))
            m_globalError = 2*m_errorThreshold;

//        if(m_myRank == 0)
//            cout << maxU << "\t" <<  avgU << "\t" <<  maxU/avgU << "\t" << m_globalError << endl;
    }
    m_counter++;
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
    MPI_Allreduce(MPI_IN_PLACE, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    m_c = 0;
    if(denominator > 0)
    {
        if(numerator/denominator > 0)
        {
            m_c = 2*sqrt(numerator/denominator);
        }
    }

    if(m_c >= 2.0)
    {
        m_c = 1.99;
    }
}
//------------------------------------------------------------------------------
void ADR::staticModifiers()
{
    for(Modifier *modifier:m_boundaryModifiers)
    {
        modifier->staticEvaluation();
    }
}
//------------------------------------------------------------------------------
void ADR::calculateStableMass()
{
#if USE_MPI
    m_mainGrid->clearGhostParticles();
    exchangeInitialGhostParticles(*m_mainGrid, *m_particles);
#endif
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
void ADR::updateGridAndCommunication()
{
    m_mainGrid->clearParticles();
    updateGrid(*m_mainGrid, *m_particles, true);

#if USE_MPI
    exchangeGhostParticles(*m_mainGrid, *m_particles);

    // Updating lists in modifiers
    int counter = 0;
    for(Modifier *modifier:m_boundaryModifiers)
    {
        updateModifierLists(*modifier, *m_particles, counter);
        counter++;
    }
#endif
}
//------------------------------------------------------------------------------
}

