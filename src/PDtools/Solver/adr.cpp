#include "adr.h"

#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/Modfiers/modifier.h"
#include "PDtools/SavePdData/savepddata.h"

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

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void ADR::stepForward(int i)
{
    save(i);
    modifiersStepOne();

    m_globalError = 2*m_errorThreshold;
    iterate(m_maxSteps);

    cout << "i = " << i << " " << m_globalError << endl;

    modifiersStepTwo();
    m_t += m_dt;
}
//------------------------------------------------------------------------------
void ADR::iterate(int maxNumberOfSteps)
{
    int counter = 0;
    int minSteps = 5;
    int nFractureRelaxationSteps = m_maxStepsFracture;
    int nParticles = m_particles->nParticles();

    do
    {
        integrateStepOne();
        zeroForcesAndStress();

        updateGridAndCommunication();
        calculateForces();
        staticModifiers();

        integrateStepTwo();

        counter++;
        if(counter > maxNumberOfSteps)
            break;
    }
    while( (m_globalError > m_errorThreshold) || (counter < minSteps));

    cout << "counter = " << counter << " error = " << m_globalError << endl;
//    calculateForces();

    //--------------------------------------------------
    if(!m_qsModifiers.empty())
    {
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        for(int i=0; i<nParticles; i++)
        {
            pair<int, int> id(i, i);
            for(Modifier *modifier:m_qsModifiers)
            {
                modifier->evaluateStepOne(id);
            }

        }
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        for(int i=0; i<nParticles; i++)
        {
            pair<int, int> id(i, i);
            for(Modifier *modifier:m_qsModifiers)
            {
                modifier->evaluateStepTwo(id);
            }

        }
        for(Modifier *modifier:m_qsModifiers)
        {
            modifier->evaluateStepTwo();
        }

        bool continueState = false;

        for(Modifier *modifier:m_qsModifiers)
        {
            if(modifier->state())
                continueState = true;
        }

        if(continueState)
        {
            m_globalError = 2*m_errorThreshold;
            iterate(nFractureRelaxationSteps);
        }
    }
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

    //-----------------------------
    calculateForces();
    calculateStableMass();
    Solver::initialize();
}
//------------------------------------------------------------------------------
void ADR::calculateStableMass()
{
    arma::vec &stableMass = m_particles->stableMass();
    int nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        pair<int, int> id(i, i);
        double stableMassMax = 0;
        for(Force *oneBodyForce:m_oneBodyForces)
        {
            double sm = oneBodyForce->calculateStableMass(id, m_dt);
            stableMassMax = max(stableMassMax, sm);
        }
        stableMass(id.second) = stableMassMax;
    }
}
//------------------------------------------------------------------------------
void ADR::integrateStepOne()
{
    using namespace arma;

    double alpha = (2.*m_dt)/(2. + m_c*m_dt);
    double beta = (2. - m_c*m_dt)/(2. + m_c*m_dt);
    double rho = sqrt(beta);
    double errorScaling = rho/(1.0001 - rho);

    mat & r = m_particles->r();
    mat & v = m_particles->v();
    mat & F = m_particles->F();
    mat & Fold = m_particles->Fold();
    arma::imat & isStatic = m_particles->isStatic();
    const vec & stableMass = m_particles->stableMass();
    double error = 0.0;
    double Fn = 0;
    double deltaFn = 0;
    int nParticles = m_particles->nParticles();
//#ifdef USE_OPENMP
//# pragma omp parallel for reduction(max: error)
//#endif
    for(int i=0; i<nParticles; i++)
    {
        int col_i = i;
        if(isStatic(col_i))
            continue;

        double drSquared = 0;

        for(int d=0; d<m_dim; d++){
            double rPrev = r(col_i, d);

            v(col_i, d) = beta*v(col_i, d) + alpha*F(col_i, d)/stableMass(col_i);
            r(col_i, d) += v(col_i, d)*m_dt;

            deltaFn += pow(Fold(col_i, d) - F(col_i, d), 2);
            Fn += pow(F(col_i, d), 2);

            Fold(col_i, d) = F(col_i, d);

            drSquared += (r(col_i, d) - rPrev)*(r(col_i, d) - rPrev);
        }

        double error_i = errorScaling*sqrt(drSquared);
        error = std::max(error, error_i);
    }

    m_globalError = error;
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
    int nParticles = m_particles->nParticles();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        int col_i = i;
        if(isStatic(col_i))
            continue;
        for(int d=0; d<m_dim; d++)
        {
            numerator += -v(col_i, d)*(F(col_i, d) - Fold(col_i, d))
                    /(stableMass(col_i)*m_dt);
            denominator += v(col_i, d)*v(col_i, d);
        }
    }

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
}

