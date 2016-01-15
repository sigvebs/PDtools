#include "simplefracture.h"

#include "PDtools/Particles/pd_particles.h"
namespace PDtools
{
//------------------------------------------------------------------------------
SimpleFracture::SimpleFracture(double alpha):
    m_alpha(4*alpha)
{

}
//------------------------------------------------------------------------------
SimpleFracture::~SimpleFracture()
{

}
//------------------------------------------------------------------------------
void SimpleFracture::initialize()
{
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexVolume = m_particles->getParamId("volume");
    m_indexMicromodulus= m_particles->getParamId("micromodulus");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_indexDr0 = m_particles->getPdParamId("dr0");
    m_indexForceScaling = m_particles->getPdParamId("forceScalingBond");
    m_idToCol = &m_particles->idToCol();
    m_data = &m_particles->data();

    m_broken = false;
    m_state = false;
}
//------------------------------------------------------------------------------
void SimpleFracture::evaluateStepOne(const int id_i, const int i)
{
    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;
    const double vol_i = (*m_data)(i, m_indexVolume);
    const double c_i = (*m_data)(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if((*m_data)(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const double c_j = (*m_data)(j, m_indexMicromodulus);
        const double g_ij = con.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        const double vol_j = (*m_data)(j, m_indexVolume);
        const double vol = vol_i + vol_j;
        const double s = con.second[m_indexStretch];
        const double dr0 = con.second[m_indexDr0];
        const double sc = sqrt(m_alpha / (c_ij*dr0*vol));

        if(s > sc)
        {
            con.second[m_indexConnected] = 0;
            m_broken = true;
        }
    }
}
//------------------------------------------------------------------------------
void SimpleFracture::evaluateStepTwo()
{
    if(m_broken)
    {
        m_state = true;
    }
    else
    {
        m_state = false;
    }

    m_broken = false;
}
//------------------------------------------------------------------------------
}
