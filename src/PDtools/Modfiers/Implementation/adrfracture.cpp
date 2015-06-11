#include "adrfracture.h"

#include "PDtools/Force/force.h"
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
void ADRfracture::initialize()
{
    m_indexS0 = m_particles->getParamId("s0");
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexS00 = m_particles->registerPdParameter("s00");
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
        }
    }

    m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
    m_maxStretch = std::numeric_limits<double>::min();
    m_state = false;
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepOne(const pair<int, int> &pIdcol)
{
    int id_i = pIdcol.first;
    int col_i = pIdcol.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    double s0_i = (*m_data)(col_i, m_indexS0);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    double s0_new = numeric_limits<double>::min();

    for(auto &con:PDconnections)
    {
        int id_j = con.first;
        int col_j = (*m_pIds)[id_j];
        if((*m_data)(col_j, m_indexUnbreakable) >= 1)
            continue;

        double s = con.second[m_indexStretch];
        double s00 = con.second[m_indexS00];
        double s0_j = (*m_data)(col_j, m_indexS0);
        double s0 = std::min(s0_i, s0_j);

        if(s > s0)
        {
            if(s > m_maxStretch)
            {
                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
                m_maxStretch = s;
            }
        }

        double s0_tmp = s00;
        if(s < 0)
        {
            s0_tmp -= m_alpha*s;
        }
        s0_new = std::max(s0_new, s0_tmp);
    }

    (*m_data)(col_i, m_indexS_tmp) = s0_new;
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepTwo(const pair<int, int> &pIdcol)
{
    int col_i = pIdcol.second;
    (*m_data)(col_i, m_indexS0) = (*m_data)(col_i, m_indexS_tmp);
}
//------------------------------------------------------------------------------
void ADRfracture::evaluateStepTwo()
{
    if(m_maxPId.first != -1)
    {
        m_state = true;
        int pId = m_maxPId.first;
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
        PDconnections.erase(remove(begin(PDconnections), end(PDconnections), *m_maxPId.second),
                             end(PDconnections) );
    }
    else
    {
        m_state = false;
    }

    m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
    m_maxStretch = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
}
