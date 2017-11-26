#include "pd_bondforce_hourglass.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondForce_hourglass::PD_bondForce_hourglass(PD_Particles &particles):
    PD_bondForce(particles)
{
//    m_iC_hg = m_particles.registerParameter("micromodulus_hg", 1);
//    m_particles.setNeedGhostR0(true);
//    _F = zeros(m_dim, m_dim);

//    switch(m_dim) {
//    case 1:
//        n_deformationGradient = 1;
//        m_indexStress[0] = m_particles.registerParameter("F_x");
//        break;
//    case 2:
//        n_deformationGradient = 3;
//        m_iF[0] = m_particles.registerParameter("F_x");
//        m_iF[1] = m_particles.registerParameter("F_y");
//        m_iF[2] = m_particles.registerParameter("F_xy");
//        break;
//    }
}
//------------------------------------------------------------------------------
void PD_bondForce_hourglass::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    PD_bondForce::initialize(E, nu, delta, dim, h, lc);

}
//------------------------------------------------------------------------------
void PD_bondForce_hourglass::calculateForces(const int id_i, const int i)
{
//    const double c_i = m_data(i, m_indexMicromodulus);
//    vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id_i);

//    _F.zeros();

//    const int nConnections = PDconnections_i.size();
//    double dr_ij[m_dim];
//    double dr0_ij[m_dim];

//    for(int l_j=0; l_j<nConnections; l_j++) {
//        auto &con_i = PDconnections_i[l_j];
//        if(con_i.second[m_indexConnected] <= 0.5)
//            continue;

//        const int id_j = con_i.first;
//        const int j = m_idToCol.at(id_j);

//        const double c_j = m_data(j, m_indexMicromodulus);
//        const double vol_j = m_data(j, m_indexVolume);
//        const double dr0 = con_i.second[m_indexDr0];
//        const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
//        const double g_ij = con_i.second[m_indexForceScaling];
//        const double c_ij = 0.5*(c_i + c_j)*g_ij;

//        double dr2 = 0;
//        for(int d=0; d<m_dim; d++) {
//            dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
//            dr_ij[d] = m_r(j, d) - m_r(i, d);
//            dr2 += dr_ij[d]*dr_ij[d];
//        }

//        const double dr = sqrt(dr2);
//        const double ds = dr - dr0;
//        const double s = ds/dr0;
//        const double fbond_ij = c_ij*s*vol_j*volumeScaling_ij/dr;

//        for(int d=0; d<m_dim; d++) {
//            m_F(i, d) += dr_ij[d]*fbond_ij;

//            for(int d2=0; d2<m_dim; d2++) {
//                _F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
//            }
//        }

//        con_i.second[m_indexStretch] = s;
//    }
}
//------------------------------------------------------------------------------
}
