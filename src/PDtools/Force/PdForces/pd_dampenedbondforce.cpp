#include "pd_dampenedbondforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_dampenedBondForce::PD_dampenedBondForce(PD_Particles &particles, double dt, double c):
    PD_bondForce(particles),
    m_dt(dt),
    m_c(c)
{

}
//------------------------------------------------------------------------------
void PD_dampenedBondForce::calculateForces(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
    const int nParticles = m_particles.nParticles();
#endif
    vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections_i.size();
    double dr_ij[m_dim];

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con_i = PDconnections_i[l_j];
        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con_i.first;
        const int j = m_idToCol.at(id_j);

#if USE_N3L
        if(j<i)
            continue;
#endif
        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con_i.second[m_indexDr0];
        const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
        const double g_ij = con_i.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;
        const double s_prev = con_i.second[m_indexStretch];

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        const double s = ds/dr0;
        const double dsdt = m_c*(s - s_prev)/m_dt;
        const double s_bond = s + dsdt;
        const double fbond_ij = c_ij*s_bond*vol_j*volumeScaling_ij/dr;

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*fbond_ij;
        }

        con_i.second[m_indexStretch] = s;
#if USE_N3L
        if(j > i && j < nParticles)
        {
            const int myPos_j = con_i.second[m_indexMyPdPosition];
            vector<pair<int, vector<double>>> & PDconnections_j = m_particles.pdConnections(id_j);
            auto &con_j = PDconnections_j[myPos_j];
            const double volumeScaling_ji = con_j.second[m_indexVolumeScaling];
            const double fbond_ji = -c_ij*s*vol_i*volumeScaling_ji/dr;
            for(int d=0; d<m_dim; d++)
            {
                m_F(j, d) += dr_ij[d]*fbond_ji;
            }
            con_j.second[m_indexStretch] = s;
        }
#endif
    }
}
//------------------------------------------------------------------------------
}
