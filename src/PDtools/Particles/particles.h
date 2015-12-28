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
    int m_dim = 0;
    mat m_r;
    unordered_map<int, int> m_pIds;
    ivec m_posToId;
    ivec m_isStatic;
    unordered_map<string, int> m_parameters;
    mat m_data;

    int m_verletUpdateFreq = 30;
    unordered_map<string , int> m_verletListIds;
    vector<unordered_map<int, vector<int>>> m_verletLists;

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
    dim() const
    {
        return m_dim;
    }

    void
    dim(int d);

    mat &
    r()
    {
        return m_r;
    }

    unordered_map<int, int> &
    pIds()
    {
        return m_pIds;
    }

    ivec &
    get_id()
    {
        return m_posToId;
    }

    mat &
    data()
    {
        return m_data;
    }

    ivec &
    isStatic()
    {
        return m_isStatic;
    }

    unordered_map<string, int> &
    parameters();
    int &
    parameters(string id);
    int
    registerVerletList(string verletId);
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
};
//------------------------------------------------------------------------------
}
#endif // PARTICLES_H
