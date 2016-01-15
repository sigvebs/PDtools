#include "micropolarfracture.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
MicropolarFracture::MicropolarFracture(double theta):
    m_thetaCritical(theta)
{

}
//------------------------------------------------------------------------------
MicropolarFracture::~MicropolarFracture()
{

}
//------------------------------------------------------------------------------
void MicropolarFracture::registerParticleParameters()
{
    m_indexS0 = m_particles->registerParameter("s0");
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_indexTheta = m_particles->registerPdParameter("theta");
    m_data = &m_particles->data();
    m_initialGhostParameters = {"s0"};
    m_ghostParameters = {"s0"};
    m_idToCol = &m_particles->idToCol();
    m_neededProperties  = {pair<string, int>("PdAngle", 1)};
}
//------------------------------------------------------------------------------
void MicropolarFracture::initialize()
{
}
//------------------------------------------------------------------------------
void MicropolarFracture::evaluateStepOne(const int id_i, const int i)
{
    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;

    const double s0_i = (*m_data)(i, m_indexS0);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if((*m_data)(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const double stretch = con.second[m_indexStretch];
        const double s0_j = (*m_data)(j, m_indexS0);
        const double s0 = std::min(s0_i, s0_j);
        const double theta = fabs(con.second[m_indexTheta]);

        cout << stretch << endl;
        if(stretch > s0)
        {
            con.second[m_indexConnected] = 0;
        }
//        if(stretch > s0 || theta > m_thetaCritical)
//        {
//            con.second[m_indexConnected] = 0;
//        }
    }
}
//------------------------------------------------------------------------------
}

