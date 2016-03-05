#ifndef PARTICLES_H
#define PARTICLES_H

#include <assert.h>
#include <unordered_map>
#include <armadillo>

namespace PDtools
{
using namespace std;
using namespace arma;
//------------------------------------------------------------------------------
// Storing all particle data
//------------------------------------------------------------------------------

class Particles
{
protected:
    string m_type;
    unsigned int m_nParticles = 0;
    unsigned int m_totParticles = 0;
    unsigned int m_maxParticles = 1e8;
    int m_dim = 0;

    // All particles
    mat m_r;
    mat m_v;
    mat m_data;
    ivec m_colToId;
    ivec m_isStatic;
    unordered_map<int, int> m_idToCol;

    // General properties
    unordered_map<string, int> m_parameters;
    int m_verletUpdateFreq = 30;
    unordered_map<string , int> m_verletListIds;
    vector<unordered_map<int, vector<int>>> m_verletLists;
    vector<string> m_ghostParametersString;
    vector<int> m_ghostParameters;
    int m_nGhostParticles = 0;
    int m_newId;
    int m_needGhostVelocity = 0;

    enum ErrorCodes
    {
        ZeroParticles,
        OutOfBounds,
        ParameterExist,
        ParameterDoesNotExist
    };
public:
    Particles();
    virtual
    ~Particles();
    virtual void
    initializeMatrices();
    const string &
    type() const;
    void
    type(string t);
    unsigned int
    nParticles() const;
    void
    nParticles(int nP);
    int
    dim() const;
    void
    maxParticles(int mp);
    void
    totParticles(int mp);
    unsigned int
    totParticles() const;
    void
    dim(int d);
    void
    nGhostParticles(int ngp);
    int
    nGhostParticles();
    mat &
    r();
    mat &
    v();
    unordered_map<int, int> &
    idToCol();

    ivec &
    colToId();

    mat &
    data();

    ivec &
    isStatic();

    virtual void
    deleteParticleById(const int deleteId);

    unordered_map<string, int> &
    parameters();

    int parameters(const string &id);

    int
    registerVerletList(const string &verletId);

    int
    getVerletSize() const;

    int
    getVerletId(string verletId) const;

    vector<int> &
    verletList(int pId, int verletId=0);

    unordered_map<int, vector<int>> &
    wholeVerletList(int verletId=0);

    void
    setVerletList(int pId, vector<int> vList, int verletId=0);
    void
    clearVerletList(int verletId);
    bool
    hasParameter(string paramId);
    int
    getParamId(string paramId);
    void
    setParameter(string paramId, double value);
    int
    registerParameter(string paramId, double value=0);
    void
    scaleParameter(const string &paramId, double value);
    int
    verletUpdateFreq() const;
    void
    setVerletUpdateFreq(int verletUpdateFreq);
    void
    addGhostParameter(const string &g_parameter);
    const vector<int> &
    ghostParameters();
    const vector<string> &
    ghostParametersString();

    int newId();
    int needGhostVelocity() const;
    void setNeedGhostVelocity(int needGhostVelocity);
};
//------------------------------------------------------------------------------
// Inline functions
inline int
Particles::dim() const
{
    return m_dim;
}

inline void
Particles::maxParticles(int mp)
{
    m_maxParticles = mp;
}

inline void
Particles::totParticles(int mp)
{
    m_totParticles = mp;
}

inline unsigned int
Particles::totParticles() const
{
    return m_totParticles;
}

inline mat &
Particles::r()
{
    return m_r;
}

inline mat &
Particles::v()
{
    return m_v;
}

inline unordered_map<int, int> &
Particles::idToCol()
{
    return m_idToCol;
}

inline ivec &
Particles::colToId()
{
    return m_colToId;
}

inline mat &
Particles::data()
{
    return m_data;
}

inline ivec &
Particles::isStatic()
{
    return m_isStatic;
}

}
#endif // PARTICLES_H
