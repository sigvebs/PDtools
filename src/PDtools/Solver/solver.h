#ifndef SOLVER_H
#define SOLVER_H

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
    vector<Modifier *> m_modifiers;
    vector<Modifier*> m_qsModifiers;
    int m_dim = 3;
    int m_steps = 0;
    double m_dt = 0;
    double m_t = 0;
    int m_saveInterval = 10;
    std::vector<std::string> m_saveParameters;
    string m_savePath = "testGeometries";
    SavePdData *saveParticles;
    double m_errorThreshold = 1.e-11;
    int m_indexStress[6];

    enum SolverErrorMessages
    {
        NumberOfStepNotSet,
        ParticlesNotSet
    };

    double m_E0 = 1.;
    double m_L0 = 1.;
    double m_v0 = 1.;
    double m_t0 = 1.;
    double m_rho0 = 1;
public:
    Solver();

    virtual ~Solver();

    virtual void
    solve() = 0;
    virtual void
    stepForward(int i) = 0;

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
    setSavePath(const string & savePath);
    void
    setSaveInterval(double saveInterval);
    void
    setSaveScaling(const double E0, const double L0, const double v0, const double t0, const double rho0);

    void
    addForce(Force* force);
    void
    addSpModifier(Modifier * modifier);
    void
    addModifier(Modifier * modifier);
    void
    addQsModifiers(Modifier * modifier);

    std::vector<std::string>
    saveParameters() const;

    void
    setSaveParameters(const std::vector<std::string> &saveParameters);

    virtual void
    updateGridAndCommunication();

    virtual void
    save(int i);

    virtual void
    initialize();

    virtual void
    modifiersStepOne();
    virtual void
    modifiersStepTwo();

    virtual void
    zeroForcesAndStress();

    void
    setADR_fracture(Modifier *ADR_fracture);
    void
    setErrorThreshold(double errorThreshold);

protected:
    void
    checkInitialization();

    virtual void
    calculateForces();

    void
    printProgress(const double progress);
};
//------------------------------------------------------------------------------
}
#endif // SOLVER_H
