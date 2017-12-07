#include "pd_lpss.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPSS::PD_LPSS(PD_Particles &particles, bool planeStress)
    : Force(particles), m_planeStress(planeStress) {
  particles.setNeedGhostR0(true);

  m_iMicromodulus = m_particles.registerParameter("micromodulus", 1);
  m_iTheta = m_particles.registerParameter("theta", 0);
  m_iThetaNew = m_particles.registerParameter("thetaNew", 0);
  m_iMass = m_particles.registerParameter("LPS_mass", 1);

  m_iVolume = m_particles.getParamId("volume");
  m_iDr0 = m_particles.getPdParamId("dr0");
  m_iVolumeScaling = m_particles.getPdParamId("volumeScaling");
  m_iConnected = m_particles.getPdParamId("connected");
  m_iBrokenNow = m_particles.registerParameter("brokenNow", 0);

  m_ghostParameters.push_back("volume");
  m_ghostParameters.push_back("theta");
  m_ghostParameters.push_back("LPS_mass");
  m_ghostParameters.push_back("micromodulus"); // For contact forces
  m_iStretch = m_particles.registerPdParameter("stretch");

  m_initialGhostParameters = {"volume", "theta", "LPS_mass"};
  m_hasUpdateState = false;
  //    m_hasUpdateState = true;

  //----------------------------------
  // Registering the rotation matrix and the shape matrix
  if (m_dim == 2) {
    m_iR[0] = m_particles.registerParameter("R_00", 1);
    m_iR[1] = m_particles.registerParameter("R_10", 0);
    m_iR[2] = m_particles.registerParameter("R_01", 0);
    m_iR[3] = m_particles.registerParameter("R_11", 1);

    m_iK[0] = m_particles.registerParameter("Ks_00");
    m_iK[1] = m_particles.registerParameter("Ks_11");
    m_iK[2] = m_particles.registerParameter("Ks_01");

    m_ghostParameters.push_back("R_00");
    m_ghostParameters.push_back("R_10");
    m_ghostParameters.push_back("R_01");
    m_ghostParameters.push_back("R_11");

    m_iStress[0] = m_particles.registerParameter("s_xx");
    m_iStress[1] = m_particles.registerParameter("s_yy");
    m_iStress[2] = m_particles.registerParameter("s_xy");
  } else if (m_dim == 3) {
    m_iR[0] = m_particles.registerParameter("R_00", 1);
    m_iR[1] = m_particles.registerParameter("R_10", 0);
    m_iR[2] = m_particles.registerParameter("R_20", 0);
    m_iR[3] = m_particles.registerParameter("R_01", 0);
    m_iR[4] = m_particles.registerParameter("R_11", 1);
    m_iR[5] = m_particles.registerParameter("R_21", 0);
    m_iR[6] = m_particles.registerParameter("R_02", 0);
    m_iR[7] = m_particles.registerParameter("R_12", 0);
    m_iR[8] = m_particles.registerParameter("R_22", 1);

    m_iK[0] = m_particles.registerParameter("Ks_00");
    m_iK[1] = m_particles.registerParameter("Ks_11");
    m_iK[2] = m_particles.registerParameter("Ks_01");
    m_iK[3] = m_particles.registerParameter("Ks_22");
    m_iK[4] = m_particles.registerParameter("Ks_02");
    m_iK[5] = m_particles.registerParameter("Ks_12");

    m_ghostParameters.push_back("R_00");
    m_ghostParameters.push_back("R_10");
    m_ghostParameters.push_back("R_20");
    m_ghostParameters.push_back("R_01");
    m_ghostParameters.push_back("R_11");
    m_ghostParameters.push_back("R_21");
    m_ghostParameters.push_back("R_02");
    m_ghostParameters.push_back("R_12");
    m_ghostParameters.push_back("R_22");

    m_iStress[0] = m_particles.registerParameter("s_xx");
    m_iStress[1] = m_particles.registerParameter("s_yy");
    m_iStress[2] = m_particles.registerParameter("s_xy");
    m_iStress[3] = m_particles.registerParameter("s_zz");
    m_iStress[4] = m_particles.registerParameter("s_xz");
    m_iStress[5] = m_particles.registerParameter("s_yz");
  }
  //----------------------------------
  m_mass = m_data.colptr(m_iMass);
  m_theta = m_data.colptr(m_iTheta);
  m_volume = m_data.colptr(m_iVolume);
  m_x = m_r.colptr(0);
  m_y = m_r.colptr(1);
  m_z = m_r.colptr(2);
  m_x0 = m_r0.colptr(0);
  m_y0 = m_r0.colptr(1);
  m_z0 = m_r0.colptr(2);

  m_Fx = m_F.colptr(0);
  m_Fy = m_F.colptr(1);
  m_Fz = m_F.colptr(2);
}
//------------------------------------------------------------------------------
void PD_LPSS::calculateForces(const int id, const int i) {
  const double theta_i = m_theta[i];
  const double m_i = m_mass[i];

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);
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
      const int j = m_idToCol_v[id_j];

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

      const double ds = sqrt(dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] +
                             dr_ij[2] * dr_ij[2]) -
                        dr0;
      con_data[m_iStretch] = ds / dr0;
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
      const int j = m_idToCol_v[id_j];

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

      // Setting the stretch
      const double ds = sqrt(dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1]) - dr0;
      con_data[m_iStretch] = ds / dr0;

      nConnected++;
    }
  }

  m_continueState = false;
}
//------------------------------------------------------------------------------
double PD_LPSS::calculatePotentialEnergyDensity(const int id_i, const int i) {
  (void)id_i;
  (void)i;

  return 0.;
}
//------------------------------------------------------------------------------
void PD_LPSS::calculatePotentialEnergy(const int id_i, const int i,
                                       int indexPotential) {
  double vol_i = m_volume[i];
  m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i) * vol_i;
}
//------------------------------------------------------------------------------
void PD_LPSS::evaluateStepOne() {
  const int nParticles = m_particles.nParticles();

// Checking for broken state
/*
for(int i=0; i<nParticles; i++) {
    const int id_i = m_colToId.at(i);
    const int broken  = m_data(i, m_iBrokenNow);
    if(broken) {
        computeMandK(id_i, i);
        m_data(i, m_iBrokenNow) = 0;
    }
}
*/
#if 1
  mat F = zeros(m_dim, m_dim);
  mat K = zeros(m_dim, m_dim);
  mat R = zeros(m_dim, m_dim);
  mat U = zeros(m_dim, m_dim);
  vec s = zeros(m_dim);
  mat V = zeros(m_dim, m_dim);

  vec dr_ij = zeros(m_dim);
  vec dr0_ij = zeros(m_dim);
  vec drr_ij = zeros(m_dim);

  for (int i = 0; i < nParticles; i++) {
    const int id = m_colToId.at(i);

    const vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(id);
    const int nConnections = PDconnections.size();

    if (m_dim == 2) {
      K(0, 0) = m_data(i, m_iK[0]);
      K(1, 1) = m_data(i, m_iK[1]);
      K(0, 1) = m_data(i, m_iK[2]);
      K(1, 0) = K(0, 1);
    }
    if (m_dim == 3) {
      K(0, 0) = m_data(i, m_iK[0]);
      K(1, 1) = m_data(i, m_iK[1]);
      K(0, 1) = m_data(i, m_iK[2]);
      K(1, 0) = K(0, 1);
      K(2, 2) = m_data(i, m_iK[3]);
      K(0, 2) = m_data(i, m_iK[4]);
      K(2, 0) = K(0, 2);
      K(1, 2) = m_data(i, m_iK[5]);
      K(2, 1) = K(1, 2);
    }

    //----------------------------------
    // Calculating the deformation gradient, F,
    // and the rotation matrix, R.
    //----------------------------------
    F.zeros();
    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      const vector<double> &con_data = con.second;

      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double vol_j = m_data(j, m_iVolume);
      const double dr0 = con.second[m_iDr0];
      const double volumeScaling_ij = con.second[m_iVolumeScaling];
      const double vol = vol_j * volumeScaling_ij;
      const double w = weightFunction(dr0);

      // Computing the deformation matrix
      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = m_r(j, d) - m_r(i, d);
        dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      }
      for (int d = 0; d < m_dim; d++) {
        for (int d2 = 0; d2 < m_dim; d2++) {
          F(d, d2) += w * dr_ij[d] * dr0_ij[d2] * vol;
        }
      }
    }
    F *= K;

    svd(U, s, V, F);
    R = U * V.t();

    if (m_dim >= 2) {
      m_data(i, m_iR[0]) = R(0, 0);
      m_data(i, m_iR[1]) = R(1, 0);
      m_data(i, m_iR[2]) = R(0, 1);
      m_data(i, m_iR[3]) = R(1, 1);
    }
    if (m_dim == 3) {
      m_data(i, m_iR[0]) = R(0, 0); // 00
      m_data(i, m_iR[1]) = R(1, 0); // 10
      m_data(i, m_iR[2]) = R(2, 0); // 20
      m_data(i, m_iR[3]) = R(0, 1); // 01
      m_data(i, m_iR[4]) = R(1, 1); // 11
      m_data(i, m_iR[5]) = R(2, 1); // 21
      m_data(i, m_iR[6]) = R(0, 2); // 02
      m_data(i, m_iR[7]) = R(1, 2); // 12
      m_data(i, m_iR[8]) = R(2, 2); // 22
    }

    //----------------------------------
    // Calculating the dilation, theta
    //----------------------------------
    const double m_i = m_mass[i];
    double theta = 0;

    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      const vector<double> &con_data = con.second;

      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double vol_j = m_data(j, m_iVolume);
      const double dr0 = con.second[m_iDr0];
      const double volumeScaling_ij = con.second[m_iVolumeScaling];
      const double vol = vol_j * volumeScaling_ij;
      const double w = weightFunction(dr0);

      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = m_r(j, d) - m_r(i, d);
        dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      }
      drr_ij = R * dr0_ij;
      double dudr0 = 0;
      for (int d = 0; d < m_dim; d++) {
        dudr0 += (dr_ij[d] - drr_ij[d]) * drr_ij[d];
      }
      theta += w * dudr0 * vol;
    }
    theta *= m_dim / m_i;
    m_theta[i] = theta;
  }
#endif
}
//------------------------------------------------------------------------------
void PD_LPSS::evaluateStepOne(const int id_i, const int i) { (void)id_i; (void) i;}
//------------------------------------------------------------------------------
double PD_LPSS::calculateStableMass(const int id_a, const int a, double dt) {
  dt *= 1.1;

  const double alpha = 2 * m_dim * m_mu;
  const double m_a = m_data(a, m_iMass);
  const arma::mat &R0 = m_particles.r0();

  double m[m_dim];
  double dr0[m_dim];
  for (int d = 0; d < m_dim; d++) {
    m[d] = 0;
  }

  const vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_a);

  double k[m_dim];

  for (int i = 0; i < m_dim; i++) {
    for (int d = 0; d < m_dim; d++) {
      k[d] = 0;
    }

    for (auto &con : PDconnections) {
      if (con.second[m_iConnected] <= 0.5)
        continue;

      const int id_b = con.first;
      const int b = m_idToCol_v[id_b];

      for (int d = 0; d < m_dim; d++) {
        dr0[d] = R0(a, d) - R0(b, d);
      }

      const double dr0Len = con.second[m_iDr0];
      const double m_b = m_data(b, m_iMass);
      const double vol_b = m_data(b, m_iVolume);
      const double volumeScaling = con.second[m_iVolumeScaling];
      const double Vb = vol_b * volumeScaling;
      const double w = weightFunction(dr0Len);

      // Check this
      const double dr0Len2 = pow(dr0Len, 2);
      double C = alpha * (1. / m_a + 1. / m_b);
      C *= w * Vb / dr0Len2;

      double sum = 0;

      for (int j = 0; j < m_dim; j++) {
        sum += fabs(dr0[j]);
      }

      k[i] += fabs(dr0[i]) * C * sum;
    }
    m[i] = k[i];
  }

  double stiffness = 0;

  for (int d = 0; d < m_dim; d++) {
    if (m[d] > stiffness) {
      stiffness = m[d];
    }
  }

  return 2. * 0.25 * pow(dt, 2) * stiffness;
}
//------------------------------------------------------------------------------
void PD_LPSS::calculateWeightedVolume() {
  const int nParticles = m_particles.nParticles();
  mat K = zeros(m_dim, m_dim);
  double dr0_ij[m_dim];

  for (int i = 0; i < nParticles; i++) {
    const int id_i = m_colToId.at(i);
    int nActiveConnections = 0;

    const vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();
    double m = 0;
    K.zeros();

    for (int l_j = 0; l_j < nConnections; l_j++) {
      auto &con = PDconnections[l_j];
      if (con.second[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];
      const double volumeScaling = con.second[m_iVolumeScaling];
      const double vol_j = m_volume[j];
      const double vol = vol_j * volumeScaling;
      const double dr0 = con.second[m_iDr0];
      const double w = weightFunction(dr0);

      m += w * dr0 * dr0 * vol;

      // Computing the shape matrix
      for (int d = 0; d < m_dim; d++) {
        dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      }
      for (int d = 0; d < m_dim; d++) {
        for (int d2 = 0; d2 < m_dim; d2++) {
          K(d, d2) += w * dr0_ij[d] * dr0_ij[d2] * vol;
        }
      }

      nActiveConnections++;
    }

    if (m_analyticalM) {
      // Remember to use the correct weightfunction
      if (m_dim == 3) {
        //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
        m = 4. / 3. * M_PI * pow(m_delta, 3); //  For w = 1/dr0^2
      } else {
        //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
        m = m_h * M_PI * pow(m_delta, 2); // For w = 1/dr0^2
      }
    }

    if (nActiveConnections <= 3) {
      K.zeros();
      cout << "WARNING--------------------------------------------------"
           << endl;
      if (m_dim == 3) {
        //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
        m = 4. / 3. * M_PI * pow(m_delta, 3); //  For w = 1/dr0^2
      } else {
        //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
        m = m_h * M_PI * pow(m_delta, 2); // For w = 1/dr0^2
      }
    } else {
      K = inv(K);
    }

    m_data(i, m_iMass) = m;

    if (m_dim == 2) {
      m_data(i, m_iK[0]) = K(0, 0);
      m_data(i, m_iK[1]) = K(1, 1);
      m_data(i, m_iK[2]) = K(0, 1);

      // Setting the initial rotation matrix
      m_data(i, m_iR[0]) = 1.; // 00
      m_data(i, m_iR[1]) = 0.; // 10
      m_data(i, m_iR[2]) = 0.; // 01
      m_data(i, m_iR[3]) = 1.; // 11
    }
    if (m_dim == 3) {
      m_data(i, m_iK[0]) = K(0, 0);
      m_data(i, m_iK[1]) = K(1, 1);
      m_data(i, m_iK[2]) = K(0, 1);
      m_data(i, m_iK[3]) = K(2, 2);
      m_data(i, m_iK[4]) = K(0, 2);
      m_data(i, m_iK[5]) = K(1, 2);

      // Setting the initial rotation matrix
      m_data(i, m_iR[0]) = 1.; // 00
      m_data(i, m_iR[1]) = 0.; // 10
      m_data(i, m_iR[2]) = 0.; // 20
      m_data(i, m_iR[3]) = 0.; // 01
      m_data(i, m_iR[4]) = 1.; // 11
      m_data(i, m_iR[5]) = 0.; // 21
      m_data(i, m_iR[6]) = 0.; // 02
      m_data(i, m_iR[7]) = 0.; // 12
      m_data(i, m_iR[8]) = 1.; // 22
    }

    // Setting the micromodulus for contact forces
    m_data(i, m_iMicromodulus) = 2 * m_dim * m_mu / m;
  }
}
//------------------------------------------------------------------------------
void PD_LPSS::computeMandK(int id_i, int i) {
  //    cout << "Reaclulcating m and K: " << id_i << ", " << i << endl;
  mat K = zeros(m_dim, m_dim);
  double dr0_ij[m_dim];

  int nActiveConnections = 0;

  const vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int nConnections = PDconnections.size();
  double m = 0;

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double vol_j = m_volume[j];
    const double vol = vol_j * volumeScaling;
    const double dr0 = con.second[m_iDr0];
    const double w = weightFunction(dr0);

    m += w * dr0 * dr0 * vol;

    // Computing the shape matrix
    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
    }
    for (int d = 0; d < m_dim; d++) {
      for (int d2 = 0; d2 < m_dim; d2++) {
        K(d, d2) += w * dr0_ij[d] * dr0_ij[d2] * vol;
      }
    }

    nActiveConnections++;
  }

  if (m_analyticalM) {
    // Remember to use the correct weightfunction
    if (m_dim == 3) {
      //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
      m = 4. / 3. * M_PI * pow(m_delta, 3); //  For w = 1/dr0^2
    } else {
      //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
      m = m_h * M_PI * pow(m_delta, 2); // For w = 1/dr0^2
    }
  }

  if (nActiveConnections <= 3) {
    K.zeros();
    cout << "WARNING--------------------------------------------------" << endl;
    if (m_dim == 3) {
      //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
      m = 4. / 3. * M_PI * pow(m_delta, 3); //  For w = 1/dr0^2
    } else {
      //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
      m = m_h * M_PI * pow(m_delta, 2); // For w = 1/dr0^2
    }
  } else {
    K = inv(K);
  }

  m_data(i, m_iMass) = m;

  if (m_dim == 2) {
    m_data(i, m_iK[0]) = K(0, 0);
    m_data(i, m_iK[1]) = K(1, 1);
    m_data(i, m_iK[2]) = K(0, 1);
  }
  if (m_dim == 3) {
    m_data(i, m_iK[0]) = K(0, 0);
    m_data(i, m_iK[1]) = K(1, 1);
    m_data(i, m_iK[2]) = K(0, 1);
    m_data(i, m_iK[3]) = K(2, 2);
    m_data(i, m_iK[4]) = K(0, 2);
    m_data(i, m_iK[5]) = K(1, 2);
  }
}
//------------------------------------------------------------------------------
void PD_LPSS::initialize(double E, double nu, double delta, int dim, double h,
                         double lc) {
  Force::initialize(E, nu, delta, dim, h, lc);
  m_delta = delta;
  m_dim = dim;
  m_nu = nu;
  m_mu = 0.5 * E / (1 + nu);
  m_k = E / (3. * (1. - 2. * nu));

  calculateWeightedVolume();
}
//------------------------------------------------------------------------------
}
