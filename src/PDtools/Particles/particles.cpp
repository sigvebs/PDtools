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
void Particles::addGhostParameter(const string &g_parameter)
{
    const int paramId = getParamId(g_parameter);
    for(const string & gp:m_ghostParametersString)
    {
        if(g_parameter == gp)
        {
            return;
        }
    }

    m_ghostParametersString.push_back(g_parameter);
    m_ghostParameters.push_back(paramId);
}
//------------------------------------------------------------------------------
const vector<int> &Particles::ghostParameters()
{
    return m_ghostParameters;
}
//------------------------------------------------------------------------------
const vector<string> &Particles::ghostParametersString()
{
    return m_ghostParametersString;
}
//------------------------------------------------------------------------------
int Particles::newId()
{
    return m_newId++;
}
//------------------------------------------------------------------------------
int Particles::needGhostVelocity() const
{
    return m_needGhostVelocity;
}
//------------------------------------------------------------------------------
void Particles::setNeedGhostVelocity(int needGhostVelocity)
{
    m_needGhostVelocity = needGhostVelocity;
}
//------------------------------------------------------------------------------
int Particles::getNeedGhostR0() const
{
    return m_needGhostR0;
}
//------------------------------------------------------------------------------
void Particles::setNeedGhostR0(int needGhostR0)
{
    m_needGhostR0 = needGhostR0;
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
    m_r = mat(PARTICLE_BUFFER*m_maxParticles, M_DIM);
    m_v  = mat(m_maxParticles*PARTICLE_BUFFER, M_DIM);
    m_data = mat(PARTICLE_BUFFER*m_maxParticles, PARAMETER_BUFFER);
    m_colToId = ivec(PARTICLE_BUFFER*m_maxParticles);
    m_isStatic = zeros<ivec>(PARTICLE_BUFFER*m_maxParticles);
    m_newId = m_maxParticles;
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
int Particles::nGhostParticles()
{
    return m_nGhostParticles;
}
//------------------------------------------------------------------------------
void Particles::nGhostParticles(int ngp)
{
    m_nGhostParticles = ngp;
}
//------------------------------------------------------------------------------
void Particles::deleteParticleById(const int deleteId)
{
    const int deleteCol = m_idToCol.at(deleteId);
    const int moveCol = m_nParticles - 1;
    const int moveId = m_colToId.at(moveCol);

    for(int d=0;d<M_DIM; d++)
    {
        m_r(deleteCol, d) = m_r(moveCol, d);
    }
    m_isStatic(deleteCol) = m_isStatic(moveCol);

    for(const auto & param:m_parameters)
    {
        const int pos = param.second;
        m_data(deleteCol, pos) = m_data(moveCol, pos);
    }

    m_colToId[deleteCol] = moveId;
    m_idToCol[moveId] = deleteCol;
    m_idToCol.erase(deleteId);
    m_nParticles--;
//    cout << "delete: " << deleteId << " col:" << deleteCol;
//    cout << "move: " << moveId << " col:" << moveCol;
}
//------------------------------------------------------------------------------
unordered_map<string, int> &Particles::parameters()
{
    return m_parameters;
}
//------------------------------------------------------------------------------
int Particles::parameters(const string &id)
{
    return m_parameters.at(id);
}
//------------------------------------------------------------------------------
int Particles::registerVerletList(const string &verletId)
{
    int id = m_verletListIds.size();
    m_verletListIds[verletId] = id;
    unordered_map<int, vector<int>> list;
    m_verletLists.push_back(list);
    return id;
}
//------------------------------------------------------------------------------
int Particles::getVerletSize() const
{
    return m_verletListIds.size();
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

    for(unsigned int p=0;p<m_nParticles; p++)
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
#ifdef DEBUG
        cerr << "WARNING: Parameter '" << paramId
             << "' already registered. Using that." << endl;
#endif
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

    for(unsigned int p=0;p<m_nParticles; p++)
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

    for(unsigned int p=0;p<m_nParticles; p++)
    {
        m_data(p, param_pos) *= value;
    }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
