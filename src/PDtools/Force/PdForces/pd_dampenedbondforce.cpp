#include "pd_dampenedbondforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_dampenedBondForce::PD_dampenedBondForce(PD_Particles &particles, double c):
    PD_bondForce(particles),
    m_c(c)
{
    particles.setNeedGhostVelocity(true);
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

        double dr2 = 0;
        double drdv = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        const double s = ds/dr0;
        const double dsdt = drdv/(dr0*dr);
        const double s_bond = s + m_c*dsdt;
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
void PD_dampenedBondForce::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{
    const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
    const int nParticles = m_particles.nParticles();
#endif

    const vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id_i);

    double dr_ij[m_dim];
    const int nConnections = PDconnections_i.size();

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        const auto &con_i = PDconnections_i[l_j];
        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con_i.first;
        const int j = m_idToCol.at(id_j);
#if USE_N3L // Already computed
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
        double drdv = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        const double s = ds/dr0;
        const double dsdt = drdv/(dr0*dr);
        const double s_bond = s + m_c*dsdt;
        const double bond_ij = c_ij*s_bond*vol_j*volumeScaling_ij/dr;

        m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
        m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
        }
#if USE_N3L
        if(j > i && j < nParticles)
        {
            const int myPos_j = con_i.second[m_indexMyPdPosition];
            const vector<pair<int, vector<double>>> & PDconnections_j = m_particles.pdConnections(id_j);
            const auto &con_j = PDconnections_j[myPos_j];
            const double volumeScaling_ji = con_j.second[m_indexVolumeScaling];
            const double bond_ji = c_ij*s*vol_i*volumeScaling_ji/dr;
            m_data(j, indexStress[0]) += 0.5*bond_ji*dr_ij[X]*dr_ij[X];
            m_data(j, indexStress[1]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Y];
            m_data(j, indexStress[2]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Y];

            if(m_dim == 3)
            {
                m_data(j, indexStress[3]) += 0.5*bond_ji*dr_ij[Z]*dr_ij[Z];
                m_data(j, indexStress[4]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Z];
                m_data(j, indexStress[5]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Z];
            }
        }
#endif
    }
}
//------------------------------------------------------------------------------
}
