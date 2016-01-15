#include "adrfracture.h"

//#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ADRfracture::ADRfracture(double alpha):
    m_alpha(alpha)
{
}
//------------------------------------------------------------------------------
ADRfracture::~ADRfracture()
{

}
//------------------------------------------------------------------------------
void ADRfracture::registerParticleParameters()
{
    m_indexS0 = m_particles->registerParameter("s0");
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexS00 = m_particles->registerPdParameter("s00");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_indexS_tmp = m_particles->registerParameter("s_tmp");
    m_idToCol = &m_particles->idToCol();
    m_data = &m_particles->data();

    m_initialGhostParameters = {"s0"};
    m_ghostParameters = {"s0"};
}
//------------------------------------------------------------------------------
void ADRfracture::initialize()
{
    const ivec & colToId = m_particles->colToId();

    // Setting the initial max stretch between two particles
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        const int id = colToId(i);
        const double s0_i = (*m_data)(i, m_indexS0);

        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id);
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];
            const double s0_j = (*m_data)(j, m_indexS0);
            con.second[m_indexS00] = 0.5*(s0_i + s0_j);
        }
    }

    m_maxPId = pair<int, int>(-1, -1);
    m_maxStretch = std::numeric_limits<double>::min();
    m_state = false;
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepOne(const int id_i, const int i)
{
    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;

    const double s0_i = (*m_data)(i, m_indexS0);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    double s0_new = numeric_limits<double>::min();
    int counter = 0;
    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if((*m_data)(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const double s = con.second[m_indexStretch];
        const double s00 = con.second[m_indexS00];
        const double s0_j = (*m_data)(j, m_indexS0);
        const double s0 = std::min(s0_i, s0_j);

        if(s > s0)
        {
            con.second[m_indexConnected] = 0;
            if(s > m_maxStretch)
            {
                m_maxPId = pair<int, int>(id_i, counter);
                m_maxStretch = s;
            }
        }

        double s0_tmp = s00;
        if(s < 0)
        {
            s0_tmp -= m_alpha*s;
        }
        s0_new = std::max(s0_new, s0_tmp);
        counter++;
    }

    (*m_data)(i, m_indexS_tmp) = s0_new;
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepTwo(const int id_i, const int i)
{
    (void) id_i;
    (*m_data)(i, m_indexS0) = (*m_data)(i, m_indexS_tmp);
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepTwo()
{
    if(m_maxPId.first != -1)
    {
        const int id_i = m_maxPId.first;
        const int remove = m_maxPId.second;
        m_state = true;
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
        PDconnections[remove].second[m_indexConnected] = 0;
    }
    else
    {
        m_state = false;
    }

    m_maxPId = std::pair<int, int>(-1, -1);
    m_maxStretch = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
}
