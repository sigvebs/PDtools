#ifndef SOLVER_H
#define SOLVER_H

#include <cstddef>
#include <vector>
#include <string>


using namespace std;

namespace PDtools
{
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
    const int m_dim = 3;
    int m_steps = 0;
    double m_dt = 0;
    double m_t = 0;
    int m_saveInterval = 10;
    std::vector<std::string> m_saveParameters;
    string m_savePath = "testGeometries";
    SavePdData *saveParticles;
    double m_errorThreshold = 1.e-11;

    enum SolverErrorMessages
    {
        NumberOfStepNotSet,
        ParticlesNotSet
    };

public:
    Solver();

    virtual ~Solver();

    virtual void solve() = 0;

    virtual void stepForward(int i) = 0;

    void setParticles(PD_Particles &_particles);

    void setDomain(Domain &_domain);

    void setMainGrid(Grid &_grid);

    void setSteps(int steps);

    void setDt(double _dt);

    void setT(double _t);

    void addForce(Force* force);

    void setSavePath(const string & savePath);

    void setSaveInterval(double saveInterval);

    void addSpModifier(Modifier * modifier);

    void addModifier(Modifier * modifier);

    std::vector<std::string> saveParameters() const;

    void setSaveParameters(const std::vector<std::string> &saveParameters);

    virtual void updateGridAndCommunication();

    virtual void save(int i);

    virtual void initialize();

    virtual void modifiersStepOne();
    virtual void modifiersStepTwo();

    void setADR_fracture(Modifier *ADR_fracture);
    void setErrorThreshold(double errorThreshold);

protected:
    void checkInitialization();
    virtual void calculateForces();
};
//------------------------------------------------------------------------------
}
#endif // SOLVER_H
