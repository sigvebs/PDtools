#include "computeaveragestretch.h"
#include "PDtools/Particles/pd_particles.h"

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
void ComputeAverageStretch::update(const int id_i, const int i)
{
    const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    double sAvg = 0;
    for(auto &con:PDconnections)
    {
        const double s = con.second[m_indexStretch];
        sAvg += s;
    }

    (*m_data)(i, m_indexAverageStretch) = sAvg / PDconnections.size();
}
//------------------------------------------------------------------------------
void ComputeAverageStretch::init(const int id_i, const int i)
{
    (void) id_i;
    (void) i;
}
//------------------------------------------------------------------------------
}
