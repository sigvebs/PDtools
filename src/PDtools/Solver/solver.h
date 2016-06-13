#ifndef SOLVER_H
#define SOLVER_H

#if USE_MPI
#include <mpi.h>
#endif

#include <cstddef>
#include <vector>
#include <string>

namespace PDtools
{
using namespace std;

class PD_Particles;
class Domain;
class Grid;
class Force;
class Modifier;
class SavePdData;
class CalculateProperty;

//------------------------------------------------------------------------------
class Solver
{
protected:
    Domain *m_domain = nullptr;
    PD_Particles *m_particles  = nullptr;
    Grid *m_mainGrid = nullptr;
    Modifier* m_ADR_fracture = nullptr;
    vector<Force *> m_oneBodyForces;
    vector<Modifier *> m_spModifiers;
    vector<Modifier *> m_boundaryModifiers;
    vector<Modifier*> m_qsModifiers;
    vector<CalculateProperty*> m_properties;

    int m_dim = 3;
    int m_steps = 0;
    double m_dt = 0;
    double m_t = 0;
    int m_saveInterval = 10;
    double m_errorThreshold = 1.e-11;

    SavePdData *m_saveParticles;

    int m_myRank = 0;
    int m_nCores = 1;

    enum SolverErrorMessages
    {
        NumberOfStepNotSet,
        ParticlesNotSet
    };
public:
    Solver();

    virtual ~Solver();

    virtual void
    solve() = 0;
    virtual void
    stepForward(int i) = 0;

    void
    applyBoundaryConditions();

    void
    setParticles(PD_Particles &_particles);
    void
    setDomain(Domain &_domain);
    void
    setMainGrid(Grid &_grid);
    void
    setSteps(int steps);
    void
    setDt(double _dt);
    void
    setT(double _t);
    void
    setDim(double _dim);

    void
    setSaveInterval(double saveInterval);

    void
    addForce(Force* force);
    void
    addSpModifier(Modifier * modifier);
    void
    addBoundaryModifier(Modifier * modifier);
    void
    addQsModifiers(Modifier * modifier);

    virtual void
    updateGridAndCommunication();

    virtual void
    updateGhosts();

    virtual void
    save(int timesStep);

    virtual void
    initialize();

    virtual void
    modifiersStepOne();
    virtual void
    modifiersStepTwo();

    virtual void
    zeroForces();

    void
    setADR_fracture(Modifier *ADR_fracture);
    void
    setErrorThreshold(double errorThreshold);

    void
    setRankAndCores(int rank, int cores);

    void
    setCalculateProperties(vector<CalculateProperty *> &calcProp);

    void
    setSaveParticles(SavePdData *saveParticles);

protected:
    void
    checkInitialization();

    virtual void
    calculateForces(int timeStep);

    void
    updateProperties(const int timeStep);

    void
    printProgress(const double progress);
};
//------------------------------------------------------------------------------
}
#endif // SOLVER_H
