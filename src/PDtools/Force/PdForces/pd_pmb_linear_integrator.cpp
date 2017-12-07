#include "pd_pmb_linear_integrator.h"

//------------------------------------------------------------------------------
namespace PDtools {
PD_PMB_LINEAR_INTEGRATOR::PD_PMB_LINEAR_INTEGRATOR(PD_Particles &particles,
                                                   double lc, double delta,
                                                   double alpha)
    : PD_bondForce(particles), m_lc(lc), m_delta(delta), m_alpha(alpha) {}
//------------------------------------------------------------------------------
void PD_PMB_LINEAR_INTEGRATOR::calculateForces(const int id, const int i) {
  const double c_i = m_data(i, m_indexMicromodulus);
  vector<pair<int, vector<double>>> &PDconnections_i =
      m_particles.pdConnections(id);
  const double lc_half = 0.5 * m_lc;
  const double lc_c = 2. / 6. * m_lc;

  const int nConnections = PDconnections_i.size();
  double dr_ij[m_dim];
  //    double dr0_ij[m_dim];
  //    double weights[m_dim];

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con_i = PDconnections_i[l_j];
    if (con_i.second[m_indexConnected] <= 0.5)
      continue;

    const int id_j = con_i.first;
    const int j = m_idToCol.at(id_j);

    const double c_j = m_data(j, m_indexMicromodulus);
    double vol_j = m_data(j, m_indexVolume);
    const double dr0 = con_i.second[m_indexDr0];
    const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
    vol_j *= volumeScaling_ij;
    //        const double g_ij = con_i.second[m_indexForceScaling];
    //        const double c_ij = 0.5*(c_i + c_j)*g_ij;
    const double c_ij = 0.5 * (c_i + c_j);

    double dr2 = 0;

    for (int d = 0; d < m_dim; d++) {
      //            dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
      //            weights[d] = 1.;
    }

    const double dr = sqrt(dr2);
    const double ds = dr - dr0;
//    const double e = dr / dr0 + 1;
    double s;

    if (m_delta > dr0 + lc_c) {
      s = ds * (1. / (dr0 - lc_half) + 2. / dr0) / 3.;
    } else if (m_delta > dr0 - lc_c) {
      s = ds * 1. / (dr0 - lc_half);
    } else {
      s = ds * (1. / (dr0 - lc_half) + 4. / dr0 + 1. / (dr0 + lc_half)) / 6.;
      //            s = ((dr*e - dr0 + lc_half)/(dr0 - lc_half)
      //                 + 4.*ds/dr0
      //                 + (dr*e - dr0 - lc_half)/(dr0 + lc_half) )/6.;
    }

    s = ds / dr0;

    //        const double s = ds/dr0;
    //        const double s = 0.5*ds*(1./(dr0 + lc_half) + 1./(dr0 - lc_half));
    //        const double s = ds*(1./(dr0 - lc_half) + 4./dr0 + 1./(dr0 +
    //        lc_half))/6.;
    const double fbond_ij = c_ij * s * vol_j / dr;

    for (int d = 0; d < m_dim; d++) {
      m_F(i, d) += dr_ij[d] * fbond_ij;
    }

    con_i.second[m_indexStretch] = s;
  }
}
//------------------------------------------------------------------------------
}
