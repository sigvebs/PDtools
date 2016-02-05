#include "pd_particles.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
vector<vector<int> > PD_Particles::getSendtParticles2() const
{
    return m_sendtParticles2;
}
//------------------------------------------------------------------------------
void PD_Particles::setSendtParticles2(const vector<vector<int> > &sendtParticles2)
{
    m_sendtParticles2 = sendtParticles2;
}
//------------------------------------------------------------------------------
vector<vector<int> > PD_Particles::getReceivedParticles2() const
{
    return m_receivedParticles2;
}
//------------------------------------------------------------------------------
void PD_Particles::setReceivedParticles2(const vector<vector<int> > &receivedParticles2)
{
    m_receivedParticles2 = receivedParticles2;
}
//------------------------------------------------------------------------------
PD_Particles::PD_Particles()
{
}
//------------------------------------------------------------------------------
PD_Particles::~PD_Particles()
{
}
//------------------------------------------------------------------------------
void PD_Particles::initializeMatrices()
{
    Particles::initializeMatrices();

    m_r0 = mat(m_maxParticles*PARTICLE_BUFFER, DIM);
    m_F  = mat(m_maxParticles*PARTICLE_BUFFER, DIM);
    m_stableMass = vec(m_maxParticles*PARTICLE_BUFFER);
    m_Fold  = mat(m_maxParticles*PARTICLE_BUFFER, DIM);
}
//------------------------------------------------------------------------------
void PD_Particles::initializeBodyForces()
{
    m_u = mat(m_nParticles, DIM); // TMP
    m_b = mat(m_nParticles, DIM); // TMP
}
//------------------------------------------------------------------------------
void PD_Particles::deleteParticleById(const int deleteId)
{
    const int deleteCol = m_idToCol.at(deleteId);
    const int moveCol = m_nParticles - 1;
    const int moveId = m_colToId.at(moveCol);
    m_idToCol.erase(deleteId);

    for(int d=0;d<DIM; d++)
    {
        m_r(deleteCol, d) = m_r(moveCol, d);
        m_r0(deleteCol, d) = m_r0(moveCol, d);
        m_v(deleteCol, d) = m_v(moveCol, d);
        m_F(deleteCol, d) = m_F(moveCol, d);
        m_Fold(deleteCol, d) = m_Fold(moveCol, d);
    }

    m_stableMass(deleteCol) = m_stableMass(moveCol);
    m_isStatic(deleteCol) = m_isStatic(moveCol);

    for(const auto & param:m_parameters)
    {
        const int pos = param.second;
        m_data(deleteCol, pos) = m_data(moveCol, pos);
    }
    m_PdConnections[deleteId].clear();
    m_PdConnections[deleteId] = m_PdConnections[moveId];
    m_colToId[deleteCol] = moveId;
    m_colToId[moveCol] = -1;
    m_idToCol[moveId] = deleteCol;
    m_nParticles--;
//    cout << " Ferdig: delete:" << deleteId << " delCol:" <<deleteCol << " moveCol:" << moveCol << endl;
}
//------------------------------------------------------------------------------
const unordered_map<string, int> &PD_Particles::PdParameters() const
{
    return m_PdParameters;
}
//------------------------------------------------------------------------------
int PD_Particles::getPdParamId(string paramId) const
{
    if(m_PdParameters.count(paramId) != 1)
    {
        cerr << "ERROR: accessing a PD_particles parameter that does not exist: "
             << paramId << endl;
        throw ParameterDoesNotExist;
    }
    return m_PdParameters.at(paramId);
}
//------------------------------------------------------------------------------
int PD_Particles::registerPdParameter(string paramId, double value)
{
    if(m_PdParameters.count(paramId) == 1)
    {
#ifdef DEBUG
        cerr << "WARNING: PD-parameter already registered: "
             << paramId << endl;
#endif
        return m_PdParameters.at(paramId);
    }
    int pos = m_PdParameters.size();
    m_PdParameters[paramId] = pos;

    // Adding the new parameter to all connections
    for(unsigned int col=0;col<m_nParticles; col++)
    {
        const int id = m_colToId(col);
        for(auto &con:m_PdConnections[id])
        {
            con.second.push_back(value);
        }
    }

    return pos;
}
//------------------------------------------------------------------------------
void PD_Particles::dimensionalScaling(const double E0, const double L0, const double v0,
                                      const double t0, const double rho0)
{
    (void) E0;
    (void) t0;

    vector<pair<int, double>> parameterScaling;

    if(hasParameter("rho"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("rho"), rho0));
    }
    if(hasParameter("mass"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("mass"), rho0));
    }
    if(hasParameter("volume"))
    {
        parameterScaling.push_back(pair<int, double>(getParamId("volume"), L0*L0*L0));
    }

    for(unsigned int i=0; i<m_nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            m_r(i, d) /= L0;
            m_r0(i, d) /= L0;
            m_v(i, d) /= v0;
        }

        for(auto pos_scaling:parameterScaling)
        {
            m_data(i, pos_scaling.first) /= pos_scaling.second;
        }
    }
}
//------------------------------------------------------------------------------
void PD_Particles::clearGhostParameters()
{
    m_ghostParameters.clear();
    m_ghostParametersString.clear();
}
//------------------------------------------------------------------------------
}
