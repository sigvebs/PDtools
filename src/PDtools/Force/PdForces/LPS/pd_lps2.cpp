#include "pd_lps2.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPS2::PD_LPS2(PDtools::PD_Particles &particles, bool planeStress,
                 bool analyticalM)
    : PD_LPS(particles, planeStress, analyticalM) {
  m_hasUpdateState = false;
}
//------------------------------------------------------------------------------
void PD_LPS2::calculateForces(const int id, const int i) {
  const double theta_i = m_theta[i];
  const double m_i = m_mass[i];

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);
  const int nConnections = PDconnections.size();
  double dr_ij[m_dim];

  if (m_dim == 3) {
    vector<double *> d_stress;
    for (int i = 0; i < 6; i++)
      d_stress.push_back(m_data.colptr(m_indexStress[i]));

    //----------------------------------
    // TMP - standard stres calc from
    for (int k = 0; k < 6; k++)
      d_stress[k][i] = 0;
    //----------------------------------
    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const double m_j = m_mass[j];
      const double theta_j = m_theta[j];
      const double vol_j = m_volume[j];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double volume = vol_j * volumeScaling;

      const double dr0 = con_data[m_iDr0];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      dr_ij[2] = m_z[j] - m_z[i];
      double dr2 =
          dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      const double dr = sqrt(dr2);
      const double ds = dr - dr0;
      double bond = m_c * (theta_i * m_i + theta_j * m_j) * dr0;
      bond += m_alpha * (m_i + m_j) * ds;
      bond *= w * volume / dr;

      m_Fx[i] += dr_ij[0] * bond;
      m_Fy[i] += dr_ij[1] * bond;
      m_Fz[i] += dr_ij[2] * bond;

      con_data[m_iStretch] = ds / dr0;
      //----------------------------------
      // TMP - standard stres calc from
      d_stress[0][i] += 0.5 * dr_ij[0] * dr_ij[0] * bond;
      d_stress[1][i] += 0.5 * dr_ij[1] * dr_ij[1] * bond;
      d_stress[2][i] += 0.5 * dr_ij[0] * dr_ij[1] * bond;
      d_stress[3][i] += 0.5 * dr_ij[2] * dr_ij[2] * bond;
      d_stress[4][i] += 0.5 * dr_ij[0] * dr_ij[2] * bond;
      d_stress[5][i] += 0.5 * dr_ij[1] * dr_ij[2] * bond;
      //----------------------------------
    }
  } else { // dim2
    //----------------------------------
    // TMP - standard stres calc from MD
    for (int k = 0; k < 3; k++)
      m_data(i, m_indexStress[k]) = 0;
    //----------------------------------
    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double m_j = m_mass[j];
      const double theta_j = m_theta[j];
      const double vol_j = m_volume[j];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double volume = vol_j * volumeScaling;

      const double dr0 = con_data[m_iDr0];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      double dr2 = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1];

      const double dr = sqrt(dr2);
      const double ds = dr - dr0;
      double bond = m_c * (theta_i * m_i + theta_j * m_j) * dr0;
      bond += m_alpha * (m_i + m_j) * ds;
      bond *= w * volume / dr;

      m_Fx[i] += dr_ij[0] * bond;
      m_Fy[i] += dr_ij[1] * bond;

      con_data[m_iStretch] = ds / dr0;
      //----------------------------------
      // TMP - standard stres calc from MD
      m_data(i, m_indexStress[0]) += 0.5 * dr_ij[0] * dr_ij[0] * bond;
      m_data(i, m_indexStress[1]) += 0.5 * dr_ij[1] * dr_ij[1] * bond;
      m_data(i, m_indexStress[2]) += 0.5 * dr_ij[0] * dr_ij[1] * bond;
    }
  }

  m_continueState = false;
}
//------------------------------------------------------------------------------
void PD_LPS2::evaluateStepOne() {
  const int nParticles = m_particles.nParticles();

  for (int i = 0; i < nParticles; i++) {
    const int id_i = m_colToId.at(i);

    const vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();

    const double m_i = m_mass[i];
    double theta = 0;
    int nConnected = 0;

    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      const vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double vol_j = m_volume[j];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double volume = vol_j * volumeScaling;
      const double dr0 = con_data[m_iDr0];
      const double w = weightFunction(dr0);
      const double s = con_data[m_iStretch];

      theta += w * s * dr0 * dr0 * volume;
      nConnected++;
    }

    if (nConnected <= 3) {
      m_theta[i] = 0;
    } else {
      theta = m_t * theta * m_i;
      m_theta[i] = theta;
    }
  }
}
//------------------------------------------------------------------------------
}
