#include "computeaveragestretch.h"
namespace PDtools
{
//------------------------------------------------------------------------------
ComputeAverageStretch::ComputeAverageStretch(PDtools::PD_Particles &particles):
    ComputeProperty(particles)
{
    m_indexStretch = m_particles.getPdParamId("stretch");
    m_indexAverageStretch = m_particles.registerParameter("average_stretch");
    m_data = &m_particles.data();
}
//------------------------------------------------------------------------------
ComputeAverageStretch::~ComputeAverageStretch()
{

}
//------------------------------------------------------------------------------
void ComputeAverageStretch::update(const pair<int, int> &pIdcol)
{
    int pId = pIdcol.first;
//    int col_i = pIdcol.second;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    double sAvg = 0;
    for(auto &con:PDconnections)
    {
        double s = con.second[m_indexStretch];
        sAvg += s;
    }

    (*m_data)(pIdcol.second, m_indexAverageStretch) = sAvg / PDconnections.size();
}
//------------------------------------------------------------------------------
void ComputeAverageStretch::init(const pair<int, int> &pIdcol)
{
    (void) pIdcol;
}
//------------------------------------------------------------------------------
}
