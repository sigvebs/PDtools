#ifndef MODIFIER_H
#define MODIFIER_H

#include <armadillo>

using namespace std;

namespace PDtools
{
// Forward declerations
class PD_Particles;

//------------------------------------------------------------------------------
class Modifier
{
protected:
    PD_Particles *m_particles;
    bool m_state = false;
    vector<int> m_localParticleIds;
    std::vector<std::string> m_initialGhostParameters;
    std::vector<std::string> m_ghostParameters;
    std::vector<pair<string, int>> m_neededProperties; // name and update frequency
    int m_dim;
    int m_myRank = 0;
    int m_nCores = 1;
public:
    Modifier();
    virtual ~Modifier();

    virtual void
    registerParticleParameters();
    virtual void
    evaluateStepOne(const int id, const int i);
    virtual void
    updateStepOne(const int id, const int i);
    virtual void
    evaluateStepTwo(const int id, const int i);
    virtual void
    evaluateStepOne();
    virtual void
    evaluateStepTwo();
    virtual void
    staticEvaluation();
    virtual void
    initialize();
    void
    setParticles(PD_Particles &particles);
    bool
    state();

    const std::vector<std::string> &
    initalGhostDependencies();

    const std::vector<std::string> &
    ghostDependencies();

    void
    addToList(int id);

    bool
    removeFromList(const int id);

    std::vector<std::pair<std::string, int> > neededProperties() const;
    int
    dim() const;
    void
    setDim(int dim);

};
//------------------------------------------------------------------------------
}
#endif // MODIFIER_H
