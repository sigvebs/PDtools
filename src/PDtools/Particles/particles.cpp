#include "particles.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
int Particles::verletUpdateFreq() const
{
    return m_verletUpdateFreq;
}
//------------------------------------------------------------------------------
void Particles::setVerletUpdateFreq(int verletUpdateFreq)
{
    m_verletUpdateFreq = verletUpdateFreq;
}
//------------------------------------------------------------------------------
Particles::Particles()
{

}
//------------------------------------------------------------------------------
Particles::~Particles()
{

}
//------------------------------------------------------------------------------
void Particles::initializeMatrices()
{
    if(m_nParticles == 0)
    {
        cerr << "nParticles can not be zero when initializeing particle matrices" << endl;
        throw ZeroParticles;
    }
    m_r = mat(PARTICLE_BUFFER*m_nParticles, DIM);
    m_data = mat(PARTICLE_BUFFER*m_nParticles, PARAMETER_BUFFER);
    m_posToId = ivec(PARTICLE_BUFFER*m_nParticles);
    m_isStatic = zeros<ivec>(PARTICLE_BUFFER*m_nParticles);
}
//------------------------------------------------------------------------------
const string &Particles::type() const
{
    return m_type;
}
//------------------------------------------------------------------------------
void Particles::type(string t)
{
    m_type = t;
}
//------------------------------------------------------------------------------
unsigned int Particles::nParticles() const
{
    return m_nParticles;
}
//------------------------------------------------------------------------------
void Particles::nParticles(int nP)
{
    m_nParticles = nP;
}
//------------------------------------------------------------------------------
void Particles::dim(int d)
{
    m_dim = d;
}
//------------------------------------------------------------------------------
unordered_map<string, int> &Particles::parameters()
{
    return m_parameters;
}
//------------------------------------------------------------------------------
int &Particles::parameters(string id)
{
    return m_parameters.at(id);
}
//------------------------------------------------------------------------------
int Particles::registerVerletList(string verletId)
{
    int id = m_verletListIds.size();
    m_verletListIds[verletId] = id;
    unordered_map<int, vector<int>> list;
    m_verletLists.push_back(list);
    return id;
}
//------------------------------------------------------------------------------
int Particles::getVerletId(string verletId) const
{
    return m_verletListIds.at(verletId);
}
//------------------------------------------------------------------------------
vector<int> &Particles::verletList(int pId, int verletId)
{
    return m_verletLists.at(verletId).at(pId);
}
//------------------------------------------------------------------------------
unordered_map<int, vector<int> > &Particles::wholeVerletList(int verletId)
{
    return m_verletLists.at(verletId);
}
//------------------------------------------------------------------------------
void Particles::setVerletList(int pId, vector<int> vList, int verletId)
{
    m_verletLists.at(verletId)[pId] = vList;
}
//------------------------------------------------------------------------------
void Particles::clearVerletList(int verletId)
{
    m_verletLists.at(verletId).clear();
}
//------------------------------------------------------------------------------
bool Particles::hasParameter(string paramId)
{
    if(m_parameters.count(paramId) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}
//------------------------------------------------------------------------------
int Particles::getParamId(string paramId)
{
    if(m_parameters.count(paramId) != 1)
    {
        cerr << "ERROR: accessing a particle parameter that does not exist: "
             << paramId << endl;
        throw ParameterDoesNotExist;
    }
    return m_parameters.at(paramId);
}

//------------------------------------------------------------------------------
void Particles::setParameter(string paramId, double value)
{
    int param_pos = m_parameters.at(paramId);

    for(unsigned int p=0;p<m_data.n_rows; p++)
    {
        m_data(p, param_pos) = value;
    }
}
//------------------------------------------------------------------------------
int Particles::registerParameter(string paramId, double value)
{
    // It the parameter exists return the original position.
    if(m_parameters.count(paramId) > 0)
    {
        cerr << "Parameter '" << paramId
             << "' already registerd. Using that." << endl;
        return m_parameters.at(paramId);
    }

    m_parameters[paramId] = m_parameters.size();
    int param_pos = m_parameters[paramId];

    if(param_pos>=PARAMETER_BUFFER)
    {
        cerr << "The number of parameters per particle is limited to"
             << PARAMETER_BUFFER << endl;
        throw OutOfBounds;
    }

    for(unsigned int p=0;p<m_data.n_rows; p++)
    {
        m_data(p, param_pos) = value;
    }

    return param_pos;
}
//------------------------------------------------------------------------------
void Particles::scaleParameter(const string &paramId, double value)
{
    m_parameters[paramId] = m_parameters.size();
    int param_pos = m_parameters[paramId];

    if(param_pos>=PARAMETER_BUFFER)
    {
        cerr << "The number of parameters per particle is limited to"
             << PARAMETER_BUFFER << endl;
        throw OutOfBounds;
    }

    for(unsigned int p=0;p<m_data.n_rows; p++)
    {
        m_data(p, param_pos) *= value;
    }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
