#include "pd_lpss_g.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPSS_G::PD_LPSS_G(PD_Particles &particles, double G0, bool planeStress)
    : PD_LPSS(particles, planeStress), m_G0(G0) {
  m_hasStepTwoModifier = true;
  m_iUnbreakable = particles.registerParameter("unbreakable");
  m_particles.registerParameter("damage");
}
//------------------------------------------------------------------------------
void PD_LPSS_G::calculateForces(const int id_i, const int i) {
  const double theta_i = m_theta[i];
  const double m_i = m_mass[i];

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int nConnections = PDconnections.size();

  double dr_ij[m_dim];
  double dr0_ij[m_dim];
  double drr_ij[m_dim];
  double drr_ji[m_dim];
  double dud_ij[m_dim];
  double dud_ji[m_dim];

  int nConnected = 0;

  double R1[5 * m_dim - 6]; // 4 for d=2 and 9 for d=3
  double R2[5 * m_dim - 6];

  if (m_dim == 3) {
    R1[0] = m_data(i, m_iR[0]); // 00
    R1[1] = m_data(i, m_iR[1]); // 10
    R1[2] = m_data(i, m_iR[2]); // 20
    R1[3] = m_data(i, m_iR[3]); // 01
    R1[4] = m_data(i, m_iR[4]); // 11
    R1[5] = m_data(i, m_iR[5]); // 21
    R1[6] = m_data(i, m_iR[6]); // 02
    R1[7] = m_data(i, m_iR[7]); // 12
    R1[8] = m_data(i, m_iR[8]); // 22

    for (int i = 0; i < 6; i++) {
      m_data(i, m_iStress[i]) = 0;
    }
    //----------------------------------
    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol.at(id_j);

      R2[0] = m_data(j, m_iR[0]); // 00
      R2[1] = m_data(j, m_iR[1]); // 10
      R2[2] = m_data(j, m_iR[2]); // 20
      R2[3] = m_data(j, m_iR[3]); // 01
      R2[4] = m_data(j, m_iR[4]); // 11
      R2[5] = m_data(j, m_iR[5]); // 21
      R2[6] = m_data(j, m_iR[6]); // 02
      R2[7] = m_data(j, m_iR[7]); // 12
      R2[8] = m_data(j, m_iR[8]); // 22

      const double m_j = m_mass[j];
      const double theta_j = m_theta[j];
      const double vol_j = m_volume[j];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double vol = vol_j * volumeScaling;
      const double dr0 = con_data[m_iDr0];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      dr_ij[2] = m_z[j] - m_z[i];
      dr0_ij[0] = m_x0[j] - m_x0[i];
      dr0_ij[1] = m_y0[j] - m_y0[i];
      dr0_ij[2] = m_z0[j] - m_z0[i];

      // drr_ij = R_i*dr0_ij;
      drr_ij[0] = R1[0] * dr0_ij[0] + R1[3] * dr0_ij[1] + R1[6] * dr0_ij[2];
      drr_ij[1] = R1[1] * dr0_ij[0] + R1[4] * dr0_ij[1] + R1[7] * dr0_ij[2];
      drr_ij[2] = R1[2] * dr0_ij[0] + R1[5] * dr0_ij[1] + R1[8] * dr0_ij[2];

      dud_ij[0] = dr_ij[0] - (1. + theta_i / m_dim) * drr_ij[0];
      dud_ij[1] = dr_ij[1] - (1. + theta_i / m_dim) * drr_ij[1];
      dud_ij[2] = dr_ij[2] - (1. + theta_i / m_dim) * drr_ij[2];

      // drr_ji = -R_j*dr0_ij;
      drr_ji[0] = -R2[0] * dr0_ij[0] - R2[3] * dr0_ij[1] - R2[6] * dr0_ij[2];
      drr_ji[1] = -R2[1] * dr0_ij[0] - R2[4] * dr0_ij[1] - R2[7] * dr0_ij[2];
      drr_ji[2] = -R2[2] * dr0_ij[0] - R2[5] * dr0_ij[1] - R2[8] * dr0_ij[2];

      dud_ji[0] = -dr_ij[0] - (1. + theta_j / m_dim) * drr_ji[0];
      dud_ji[1] = -dr_ij[1] - (1. + theta_j / m_dim) * drr_ji[1];
      dud_ji[2] = -dr_ij[2] - (1. + theta_j / m_dim) * drr_ji[2];

      m_Fx[i] += vol * m_dim * w *
                 ((m_k * theta_i * drr_ij[0] + 2 * m_mu * dud_ij[0]) / m_i -
                  (m_k * theta_j * drr_ji[0] + 2 * m_mu * dud_ji[0]) / m_j);
      m_Fy[i] += vol * m_dim * w *
                 ((m_k * theta_i * drr_ij[1] + 2 * m_mu * dud_ij[1]) / m_i -
                  (m_k * theta_j * drr_ji[1] + 2 * m_mu * dud_ji[1]) / m_j);
      m_Fz[i] += vol * m_dim * w *
                 ((m_k * theta_i * drr_ij[2] + 2 * m_mu * dud_ij[2]) / m_i -
                  (m_k * theta_j * drr_ji[2] + 2 * m_mu * dud_ji[2]) / m_j);

      // Setting the devaitoric bond strain
      const double e_i = (dud_ij[0] * dud_ij[0] + dud_ij[1] * dud_ij[1] +
                          dud_ij[2] * dud_ij[2]);
      const double e_j = (dud_ji[0] * dud_ji[0] + dud_ji[1] * dud_ji[1] +
                          dud_ji[2] * dud_ji[2]);
      con.second[m_iStretch] = std::max(e_i, e_j) / (dr0 * dr0);
      nConnected++;
    }
  }
  if (m_dim == 2) {             // dim2
    R1[0] = m_data(i, m_iR[0]); // 00
    R1[1] = m_data(i, m_iR[1]); // 10
    R1[2] = m_data(i, m_iR[2]); // 01
    R1[3] = m_data(i, m_iR[3]); // 11

    for (int i = 0; i < 3; i++)
      m_data(i, m_iStress[i]) = 0;

    //----------------------------------
    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol.at(id_j);

      R2[0] = m_data(j, m_iR[0]); // 00
      R2[1] = m_data(j, m_iR[1]); // 10
      R2[2] = m_data(j, m_iR[2]); // 01
      R2[3] = m_data(j, m_iR[3]); // 11

      const double m_j = m_mass[j];
      const double theta_j = m_theta[j];
      const double vol_j = m_volume[j];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double vol = vol_j * volumeScaling;
      const double dr0 = con_data[m_iDr0];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      dr0_ij[0] = m_x0[j] - m_x0[i];
      dr0_ij[1] = m_y0[j] - m_y0[i];

      // drr_ij = R_i*dr0_ij;
      drr_ij[0] = R1[0] * dr0_ij[0] + R1[2] * dr0_ij[1];
      drr_ij[1] = R1[1] * dr0_ij[0] + R1[3] * dr0_ij[1];
      dud_ij[0] = dr_ij[0] - (1. + theta_i / m_dim) * drr_ij[0];
      dud_ij[1] = dr_ij[1] - (1. + theta_i / m_dim) * drr_ij[1];
      // dud_ij = dr_ij - (1. + theta_i/m_dim)*drr_ij;

      // drr_ji = -R_j*dr0_ij;
      drr_ji[0] = -R2[0] * dr0_ij[0] - R2[2] * dr0_ij[1];
      drr_ji[1] = -R2[1] * dr0_ij[0] - R2[3] * dr0_ij[1];
      dud_ji[0] = -dr_ij[0] - (1. + theta_j / m_dim) * drr_ji[0];
      dud_ji[1] = -dr_ij[1] - (1. + theta_j / m_dim) * drr_ji[1];
      // dud_ji = -dr_ij - (1. + theta_j/m_dim)*drr_ji;

      m_Fx[i] += vol * m_dim * w *
                 ((m_k * theta_i * drr_ij[0] + 2 * m_mu * dud_ij[0]) / m_i -
                  (m_k * theta_j * drr_ji[0] + 2 * m_mu * dud_ji[0]) / m_j);
      m_Fy[i] += vol * m_dim * w *
                 ((m_k * theta_i * drr_ij[1] + 2 * m_mu * dud_ij[1]) / m_i -
                  (m_k * theta_j * drr_ji[1] + 2 * m_mu * dud_ji[1]) / m_j);

      // Setting the devaitoric bond strain
      const double e_i = (dud_ij[0] * dud_ij[0] + dud_ij[1] * dud_ij[1]);
      const double e_j = (dud_ji[0] * dud_ji[0] + dud_ji[1] * dud_ji[1]);
      con.second[m_iStretch] = std::max(e_i, e_j) / (dr0 * dr0);
      nConnected++;
    }
  }

  m_continueState = false;
  m_data(i, m_iBrokenNow) = 0;
}
//------------------------------------------------------------------------------
void PD_LPSS_G::evaluateStepTwo(const int id_i, const int i) {
  if (m_data(i, m_iUnbreakable) >= 1)
    return;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  bool broken = false;
  const double theta_i = m_theta[i];

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = m_idToCol[id_j];

    if (m_data(j, m_iUnbreakable) >= 1)
      continue;

    if (con.second[m_iConnected] <= 0.5)
      continue;

    const double e_d = con.second[m_iStretch];
    const double theta_j = m_theta[j];

    if (e_d > m_e_max) {
      if (theta_i > 0 || theta_j > 0) {
        m_data(i, m_iBrokenNow) = 1;
        con.second[m_iConnected] = 0;
        m_continueState = true;
        broken = true;
      }
    }
  }

  if (broken) {
    computeMandK(id_i, i);
  }
}
//------------------------------------------------------------------------------
void PD_LPSS_G::initialize(double E, double nu, double delta, int dim, double h,
                           double lc) {
  PD_LPSS::initialize(E, nu, delta, dim, h, lc);

  // Setting the maxiamal deviatoric strain on a bond
  if (m_dim == 3) {
    m_e_max = 4. / 3. * sqrt(m_G0 / (m_mu * m_delta));
  }
  if (m_dim == 2) {
    m_e_max = sqrt(3 * M_PI * m_G0 / (2 * m_mu * m_delta));
  }
  cout << "Stretch factor:" << m_e_max << endl;
  // Squaring here means that we can use e_d^2 instead of e_d.
  m_e_max *= m_e_max;
}
//------------------------------------------------------------------------------
}
