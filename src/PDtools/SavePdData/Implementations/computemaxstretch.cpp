#include "computemaxstretch.h"


namespace PDtools
{
//------------------------------------------------------------------------------
ComputeMaxStretch::ComputeMaxStretch(PD_Particles &particles):
    ComputeProperty(particles)
{
    m_indexStretch = m_particles.getPdParamId("stretch");
    m_indexMaxStretch = m_particles.registerParameter("max_stretch");
    m_data = &m_particles.data();
}
//------------------------------------------------------------------------------
ComputeMaxStretch::~ComputeMaxStretch()
{

}
//------------------------------------------------------------------------------
void ComputeMaxStretch::update(const pair<int, int> &pIdcol)
{
    int pId = pIdcol.first;
//    int col_i = pIdcol.second;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    double sMax = 0;
    for(auto &con:PDconnections)
    {
        double s = con.second[m_indexStretch];
        sMax = max(fabs(s), sMax);
    }
    (*m_data)(pIdcol.second, m_indexMaxStretch) = sMax;
}
//------------------------------------------------------------------------------
void ComputeMaxStretch::init(const pair<int, int> &pIdcol)
{
    (void) pIdcol;
}
//------------------------------------------------------------------------------
}

