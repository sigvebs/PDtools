#ifndef PARTICLES_H
#define PARTICLES_H

#include <armadillo>
#include <assert.h>
#include <unordered_map>

using namespace std;

namespace PDtools
{
//------------------------------------------------------------------------------
// Storing all particle data
//------------------------------------------------------------------------------

class Particles
{
protected:
    std::string m_type;
    std::size_t m_nParticles = 0;
    std::size_t m_dim = 0;
    arma::mat m_r;
    std::unordered_map<int, int> m_pIds;
    arma::ivec m_posToId;
    arma::ivec m_isStatic;
    std::unordered_map<std::string, int> m_parameters;
    arma::mat m_data;

    int m_verletUpdateFreq = 30;
    std::unordered_map<std::string , int> m_verletListIds;
    std::vector<std::unordered_map<int, std::vector<int>>> m_verletLists;

    enum ErrorCodes
    {
        ZeroParticles,
        OutOfBounds,
        ParameterExist,
        ParameterDoesNotExist
    };

public:
    Particles();

    virtual ~Particles();

    virtual void initializeMatrices();

    const string & type() const;

    void type(string t);

    size_t nParticles() const;

    void nParticles(int nP);

    size_t dim() const
    {
        return m_dim;
    }

    void dim(int d);

    arma::mat & r()
    {
        return m_r;
    }

    std::unordered_map<int, int> & pIds()
    {
        return m_pIds;
    }

    arma::ivec & get_id()
    {
        return m_posToId;
    }

    arma::mat & data()
    {
        return m_data;
    }

    arma::ivec & isStatic()
    {
        return m_isStatic;
    }

    std::unordered_map<std::string, int> & parameters();

    int & parameters(string id);

    int registerVerletList(std::string verletId);

    int getVerletId(std::string verletId) const;

    vector<int> & verletList(int pId, int verletId=0);

    std::unordered_map<int, std::vector<int>> & wholeVerletList(int verletId=0);

    void  setVerletList(int pId, vector<int> vList, int verletId=0);

    void clearVerletList(int verletId);

    bool hasParameter(string paramId);

    int getParamId(string paramId);

    void setParameter(string paramId, double value);

    int registerParameter(string paramId, double value=0);
    void scaleParameter(const std::string &paramId, double value);
    int verletUpdateFreq() const;
    void setVerletUpdateFreq(int verletUpdateFreq);
};
//------------------------------------------------------------------------------
}
#endif // PARTICLES_H
