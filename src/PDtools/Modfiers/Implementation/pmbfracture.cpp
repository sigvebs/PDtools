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
}
//------------------------------------------------------------------------------
void PmbFracture::evaluateStepOne(const pair<int, int> &id_col)
{
    int pId = id_col.first;
    int col_i = id_col.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    double s0_i = (*m_data)(col_i, m_indexS0);
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
    vector<pair<int, vector<double>> *> removeParticles;

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
            removeParticles.push_back(&con);
        }

        double s0_tmp = s00;
        if(s < 0)
        {
            s0_tmp -= m_alpha*s;
        }
        s0_new = std::max(s0_new, s0_tmp);
    }

    (*m_data)(col_i, m_indexS_tmp) = s0_new;

    for(auto &removeParticle:removeParticles)
    {
        PDconnections.erase( remove(begin(PDconnections), end(PDconnections), *removeParticle),
                             end(PDconnections) );
    }
}
//------------------------------------------------------------------------------
void PmbFracture::evaluateStepTwo(const pair<int, int> &id_col)
{
    int col_i = id_col.second;
    (*m_data)(col_i, m_indexS0) = (*m_data)(col_i, m_indexS_tmp);
}
//------------------------------------------------------------------------------
}

