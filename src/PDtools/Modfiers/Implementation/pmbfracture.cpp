#include "pmbfracture.h"

#include <limits>
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PmbFracture::PmbFracture(double alpha):
    m_alpha(alpha)
{
}
//------------------------------------------------------------------------------
PmbFracture::~PmbFracture()
{

}
//------------------------------------------------------------------------------
void PmbFracture::initialize()
{
    m_indexS0 = m_particles->getParamId("s0");
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexS00 = m_particles->registerPdParameter("s00");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_indexS_tmp = m_particles->registerParameter("s_tmp");
    m_pIds = &m_particles->pIds();
    m_data = &m_particles->data();

    // Setting the initial max stretch between two particles
    for(int i=0;i<m_particles->nParticles();i++)
    {
        int pId = i;
        int col_i = i;
        double s0_i = (*m_data)(col_i, m_indexS0);

        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = (*m_pIds)[id_j];
            double s0_j = (*m_data)(col_j, m_indexS0);
            con.second[m_indexS00] = 0.5*(s0_i + s0_j);
            m_s00 = s0_i; // TMP FIX
        }
    }
}
//------------------------------------------------------------------------------
void PmbFracture::evaluateStepOne(const pair<int, int> &id_col)
{
    const int pId = id_col.first;
    const int col_i = id_col.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    const double s0_i = (*m_data)(col_i, m_indexS0);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);

    double s0_new = numeric_limits<double>::min();

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_pIds)[id_j];

        if((*m_data)(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const double stretch = con.second[m_indexStretch];
        const double s00 = con.second[m_indexS00];
        const double s0_j = (*m_data)(j, m_indexS0);
        const double s0 = std::min(s0_i, s0_j);
//        s0 *= g_ij;

        if(stretch > s0)
        {
            con.second[m_indexConnected] = 0;
        }

        double s0_tmp = s00 - m_alpha*stretch;
        s0_new = std::max(s0_new, s0_tmp);
    }

//    if(s0_new < m_s00)
//        (*m_data)(col_i, m_indexS_tmp) = m_s00;
//    else
        (*m_data)(col_i, m_indexS_tmp) = s0_new;
}
//------------------------------------------------------------------------------
void PmbFracture::evaluateStepTwo(const pair<int, int> &id_col)
{
    int col_i = id_col.second;
    (*m_data)(col_i, m_indexS0) = (*m_data)(col_i, m_indexS_tmp);
}
//------------------------------------------------------------------------------
}

