#include "adrfractureaverage.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ADRfractureAverage::ADRfractureAverage(double alpha):
    m_alpha(alpha)
{

}
//------------------------------------------------------------------------------
ADRfractureAverage::~ADRfractureAverage()
{

}
//------------------------------------------------------------------------------
void ADRfractureAverage::initialize()
{
    m_indexS0 = m_particles->getParamId("s0");
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexS00 = m_particles->registerPdParameter("s00");
    m_indexS_avg = m_particles->registerParameter("s_avg");
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
void ADRfractureAverage::evaluateStepOne(const pair<int, int> &pIdcol)
{
    int pId = pIdcol.first;
    int col_i = pIdcol.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    double s0_i = (*m_data)(col_i, m_indexS0);
    double s_i = (*m_data)(col_i, m_indexS_avg);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);


    for(auto &con:PDconnections)
    {
        int id_j = con.first;
        int col_j = (*m_pIds)[id_j];
        if((*m_data)(col_j, m_indexUnbreakable) >= 1)
            continue;

        double s0_j = (*m_data)(col_j, m_indexS0);
        double s_j = (*m_data)(col_j, m_indexS_avg);
        double s = 0.5*(s_i + s_j);
        double s0 = std::min(s0_i, s0_j);

        if(s > s0)
        {
            if(s > m_maxStretch)
            {
                m_maxPId = pair<int, pair<int, vector<double>> *>(pId, &con);
                m_maxStretch = s;
            }
        }
    }

}
//------------------------------------------------------------------------------
void ADRfractureAverage::evaluateStepTwo(const pair<int, int> &pIdcol)
{
    int pId = pIdcol.first;
    int col_i = pIdcol.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);

    double s_avg = 0;
    for(auto &con:PDconnections)
    {
        s_avg += con.second[m_indexStretch];
    }

    (*m_data)(col_i, m_indexS_avg) = s_avg / PDconnections.size();
}
//------------------------------------------------------------------------------
void ADRfractureAverage::evaluateStepTwo()
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
