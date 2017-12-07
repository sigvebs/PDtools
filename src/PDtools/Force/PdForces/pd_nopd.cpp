#include "pd_nopd.h"

namespace PDtools {
//------------------------------------------------------------------------------

PD_NOPD::PD_NOPD(PD_Particles &particles, double phi, double C, double T,
                 bool planeStress)
    : Force(particles), m_T(T), m_planeStress(planeStress) {
  m_particles.setNeedGhostR0(true);
  m_hasStepTwoModifier = true;

  m_phi = phi * M_PI / 180.;
  m_d = tan(m_phi);
  m_C = 0.5 * C * (1. / (sqrt(m_d * m_d + 1.) + m_d));

  m_cos_theta = cos(M_PI / 2. + m_phi);
  m_sin_theta = sin(M_PI / 2. + m_phi);

  m_greenStrain = true;
  m_hasUpdateState = true;
  m_delta = 1.;
}
//------------------------------------------------------------------------------
void PD_NOPD::initialize(double E, double nu, double delta, int dim, double h,
                         double lc) {
  Force::initialize(E, nu, delta, dim, h, lc);

  m_dampCoeff = 5.e13 * E / (lc * lc);
  m_dampCoeff = 5.e11 * E / (lc * lc);

  if (m_dim == 3) {
    const double nu_t = 1. / 4.;
    const double k = m_E / (3. * (1. - 2. * nu_t));
    m_C_PMB = 18. * k / (M_PI * pow(m_delta, 4));
  } else if (m_dim == 2) {
    const double nu_t = 1. / 3.;
    const double k = m_E / (2. * (1. - nu_t));
    m_C_PMB = 12 * k / (m_h * M_PI * pow(m_delta, 3));
  } else if (m_dim == 1) {
    m_C_PMB = 2 * m_E / (m_h * m_h * pow(m_delta, 2));
  }
  m_C_hg = 0.5;
  cout << "m_C_PMB: " << m_C_PMB << endl;

  //    m_F = zeros(m_dim, m_dim);
  m_strain = zeros(m_dim, m_dim);
  m_K_i = zeros(m_dim, m_dim);
  m_PK_i = zeros(m_dim, m_dim);
  m_K_j = zeros(m_dim, m_dim);
  m_PK_j = zeros(m_dim, m_dim);
  m_DGT = zeros(m_dim, m_dim);
  f_ij = {0, 0, 0};

  m_iConnected = m_particles.getPdParamId("connected");
  m_iUnbreakable = m_particles.registerParameter("unbreakable");
  m_iVolume = m_particles.getParamId("volume");
  m_iVolumeScaling = m_particles.getPdParamId("volumeScaling");
  m_iDr0 = m_particles.getPdParamId("dr0");
  m_iBrokenNow = m_particles.registerParameter("brokenNow", 0);
  m_particles.registerParameter("damage");

  m_ghostParameters.push_back("volume");
  m_initialGhostParameters.push_back("volume");

  switch (m_dim) {
  case 1:
    m_nStressStrainElements = 1;
    m_ghostParameters.push_back("s_xx");
    m_ghostParameters.push_back("K_xx");
    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStrain[0] = m_particles.registerParameter("e_xx");
    m_indexK[0] = m_particles.registerParameter("K_xx");
    break;
  case 2:
    m_nStressStrainElements = 3;
    m_indexK[0] = m_particles.registerParameter("K_xx");
    m_indexK[1] = m_particles.registerParameter("K_yy");
    m_indexK[2] = m_particles.registerParameter("K_xy");
    m_indexK[3] = m_particles.registerParameter("K_yx");
    m_indexPK[0] = m_particles.registerParameter("PK_xx");
    m_indexPK[1] = m_particles.registerParameter("PK_yy");
    m_indexPK[2] = m_particles.registerParameter("PK_xy");
    m_indexPK[3] = m_particles.registerParameter("PK_yx");
    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStress[1] = m_particles.registerParameter("s_yy");
    m_indexStress[2] = m_particles.registerParameter("s_xy");
    m_indexStress[3] = m_particles.registerParameter("s_yx");
    m_indexF[0] = m_particles.registerParameter("F_xx");
    m_indexF[1] = m_particles.registerParameter("F_yy");
    m_indexF[2] = m_particles.registerParameter("F_xy");
    m_indexF[3] = m_particles.registerParameter("F_yx");

    m_ghostParameters.push_back("PK_xx");
    m_ghostParameters.push_back("PK_yy");
    m_ghostParameters.push_back("PK_xy");
    m_ghostParameters.push_back("PK_yx");
    break;
  case 3:
    m_nStressStrainElements = 6;
    m_ghostParameters.push_back("s_xx");
    m_ghostParameters.push_back("s_yy");
    m_ghostParameters.push_back("s_xy");
    m_ghostParameters.push_back("s_xy");
    m_ghostParameters.push_back("s_xz");
    m_ghostParameters.push_back("s_yz");
    m_ghostParameters.push_back("K_xx");
    m_ghostParameters.push_back("K_yy");
    m_ghostParameters.push_back("K_xy");
    m_ghostParameters.push_back("K_xy");
    m_ghostParameters.push_back("K_xz");
    m_ghostParameters.push_back("K_yz");
    m_indexK[0] = m_particles.registerParameter("K_xx");
    m_indexK[1] = m_particles.registerParameter("K_yy");
    m_indexK[2] = m_particles.registerParameter("K_xy");
    m_indexK[3] = m_particles.registerParameter("K_zz");
    m_indexK[4] = m_particles.registerParameter("K_xz");
    m_indexK[5] = m_particles.registerParameter("K_yz");
    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStress[1] = m_particles.registerParameter("s_yy");
    m_indexStress[2] = m_particles.registerParameter("s_xy");
    m_indexStress[3] = m_particles.registerParameter("s_zz");
    m_indexStress[4] = m_particles.registerParameter("s_xz");
    m_indexStress[5] = m_particles.registerParameter("s_yz");
    m_indexStrain[0] = m_particles.registerParameter("e_xx");
    m_indexStrain[1] = m_particles.registerParameter("e_yy");
    m_indexStrain[2] = m_particles.registerParameter("e_xy");
    m_indexStrain[3] = m_particles.registerParameter("e_zz");
    m_indexStrain[4] = m_particles.registerParameter("e_xz");
    m_indexStrain[5] = m_particles.registerParameter("e_yz");
    break;
  }

  // Computing the shape tensor
  const ivec &colToId = m_particles.colToId();
  const int nParticles = m_particles.nParticles();

  for (int i = 0; i < nParticles; i++) {
    const int id_i = colToId(i);
    computeK(id_i, i);
  }
}
//------------------------------------------------------------------------------
void PD_NOPD::updateState(int id, int i) {
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);

  const int nConnections = PDconnections.size();
  double dr0_ij[m_dim];
  double dr_ij[m_dim];

  mat m_DGT = zeros(m_dim, m_dim);
  mat K = zeros(m_dim, m_dim);
  mat E = zeros(m_dim, m_dim);
  mat P = zeros(m_dim, m_dim);
  mat S_cauchy = zeros(m_dim, m_dim);

  int nConnected = 0;

  // Calculating the deformation gradient tensor F
  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];

    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol.at(id_j);

    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double vol = vol_j * volumeScaling;
    const double w = weightFunction(dr0);

    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      dr_ij[d] = m_r(j, d) - m_r(i, d);
    }

    for (int d = 0; d < m_dim; d++) {
      for (int d2 = 0; d2 < m_dim; d2++) {
        m_DGT(d, d2) += w * dr_ij[d] * dr0_ij[d2] * vol;
      }
    }
    nConnected++;
  }

  if (m_data(i, m_iBrokenNow)) {
    computeK(id, i);
    m_data(i, m_iBrokenNow) = 0;
  }

  K(0, 0) = m_data(i, m_indexK[0]);
  K(1, 1) = m_data(i, m_indexK[1]);
  K(0, 1) = m_data(i, m_indexK[2]);
  K(1, 0) = m_data(i, m_indexK[3]);

  if (nConnected <= 3) {
    m_data(i, m_indexStrain[0]) = 0;
    m_data(i, m_indexStrain[1]) = 0;
    m_data(i, m_indexStrain[2]) = 0;
    m_data(i, m_indexStress[0]) = 0;
    m_data(i, m_indexStress[1]) = 0;
    m_data(i, m_indexStress[2]) = 0;
    return;
  } else {
    m_DGT = m_DGT * K; // K = inv(K);
    E = 0.5 * m_DGT * m_DGT.t();
    E(0, 0) -= 0.5;
    E(1, 1) -= 0.5;
  }

  // Storing the deformation gradient tensor
  m_data(i, m_indexF[0]) = m_DGT(0, 0);
  m_data(i, m_indexF[1]) = m_DGT(1, 1);
  m_data(i, m_indexF[2]) = m_DGT(0, 1);
  m_data(i, m_indexF[3]) = m_DGT(1, 0);

  // Assuming linear elasticity
  if (m_dim == 2) {
    // Constituent model, linear elastic
    // Computing the PK2 stress

    if (m_planeStress) {
      const double a = m_E / (1. - m_nu * m_nu);
      P(0, 0) = a * (E(0, 0) + m_nu * E(1, 1));
      P(1, 1) = a * (E(1, 1) + m_nu * E(0, 0));
      P(0, 1) = a * (1 - m_nu) * E(0, 1);
      P(1, 0) = P(0, 1);
    } else { // Plane strain
      const double a = m_E / ((1. + m_nu) * (1. - 2 * m_nu));
      P(0, 0) = a * ((1. - m_nu) * E(0, 0) + m_nu * E(1, 1));
      P(1, 1) = a * ((1. - m_nu) * E(1, 1) + m_nu * E(0, 0));
      P(0, 1) = a * 0.5 * (1. - 2. * m_nu) * E(0, 1);
      P(1, 0) = P(0, 1);
    }

    P = m_DGT * P;
    double J = det(m_DGT);
    S_cauchy = 1. / J * P * m_DGT.t();
    P = P * K; // K = inv(K)

    m_data(i, m_indexPK[0]) = P(0, 0);
    m_data(i, m_indexPK[1]) = P(1, 1);
    m_data(i, m_indexPK[2]) = P(0, 1);
    m_data(i, m_indexPK[3]) = P(1, 0);

    m_data(i, m_indexStress[0]) = S_cauchy(0, 0);
    m_data(i, m_indexStress[1]) = S_cauchy(1, 1);
    m_data(i, m_indexStress[2]) = S_cauchy(0, 1);
    m_data(i, m_indexStress[3]) = S_cauchy(1, 0);
  }
}
//------------------------------------------------------------------------------
void PD_NOPD::calculateForces(const int id, const int i) {
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);
  const int nConnections = PDconnections.size();
  unordered_map<int, int> &m_idToCol = m_particles.idToCol();

  m_PK_i(0, 0) = m_data(i, m_indexPK[0]);
  m_PK_i(1, 1) = m_data(i, m_indexPK[1]);
  m_PK_i(0, 1) = m_data(i, m_indexPK[2]);
  m_PK_i(1, 0) = m_data(i, m_indexPK[3]);

  m_DGT(0, 0) = m_data(i, m_indexF[0]);
  m_DGT(1, 1) = m_data(i, m_indexF[1]);
  m_DGT(0, 1) = m_data(i, m_indexF[2]);
  m_DGT(1, 0) = m_data(i, m_indexF[3]);

  vec f = zeros(m_dim);
  vec dr_ij = zeros(m_dim);
  vec dr0_ij = zeros(m_dim);

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];

    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol.at(id_j);
    const double volum_j = m_data(j, m_iVolume);
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double vol_j = volum_j * volumeScaling;

    m_PK_j(0, 0) = m_data(j, m_indexPK[0]);
    m_PK_j(1, 1) = m_data(j, m_indexPK[1]);
    m_PK_j(0, 1) = m_data(j, m_indexPK[2]);
    m_PK_j(1, 0) = m_data(j, m_indexPK[3]);

    const double dr0 = con.second[m_iDr0];
    const double w = weightFunction(dr0);
    double dr = 0;

    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr += dr_ij[d] * dr_ij[d];
    }
    dr = sqrt(dr);

    f = w * (m_PK_i + m_PK_j) * dr0_ij * vol_j;

    // Adding the hourglass model
    const double h_proj = arma::dot((-dr_ij + m_DGT * dr0_ij), dr0_ij);
    f -= m_C_hg * m_C_PMB * h_proj / (dr0 * dr0) * dr_ij * vol_j;

    for (int d = 0; d < m_dim; d++) {
      m_F(i, d) += f(d);
    }
  }
}
//------------------------------------------------------------------------------
void PD_NOPD::evaluateStepTwo(int id, int i) {
  (void)id;
  (void)i;
  /*
  if(m_data(i, m_iUnbreakable) >= 1)
      return;

  vector<pair<int, vector<double>>> & PDconnections =
  m_particles.pdConnections(id);

  if(m_dim == 2)
  {
      for(auto &con:PDconnections)
      {
          const int id_j = con.first;
          const int j = m_idToCol[id_j];

          if(m_data(j, m_iUnbreakable) >= 1)
              continue;

          if(con.second[m_iConnected] <= 0.5)
              continue;

          const double sx = 0.5*(m_data(i, m_indexStress[0]) + m_data(j,
  m_indexStress[0]));
          const double sy = 0.5*(m_data(i, m_indexStress[1]) + m_data(j,
  m_indexStress[1]));
          const double sxy = 0.5*(m_data(i, m_indexStress[2]) + m_data(j,
  m_indexStress[2]));

          const double first = 0.5*(sx + sy);
          const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

          double p_2 = first + second; // max
          double p_1 = first - second; // min

          const double shear = fabs(0.5*(p_1 - p_2)*m_sin_theta);
          const double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*m_cos_theta;

          const double criticalShear = fabs(shear) - fabs(m_C - m_d*normal);
          const double criticalTensile = p_2 - m_T;

          if(criticalShear >= 0  && normal < 0)
          {
              m_data(i, m_iBrokenNow) = 1;
              con.second[m_iConnected] = 0;
          }
          else if(criticalTensile >= 0)
          {
              m_data(i, m_iBrokenNow) = 1;
              con.second[m_iConnected] = 0;
          }
      }
  }
  */
}
//------------------------------------------------------------------------------
void PD_NOPD::computeK(int id, int i) {
  mat K = zeros(m_dim, m_dim);
  vector<pair<int, vector<double>>> &PDconnections_i =
      m_particles.pdConnections(id);
  const int nConnections = PDconnections_i.size();
  double dr0_ij[m_dim];

  int nConnected = 0;
  for (int l_j = 0; l_j < nConnections; l_j++) {

    const auto &con_i = PDconnections_i[l_j];
    if (con_i.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con_i.first;
    const int j = m_idToCol.at(id_j);

    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con_i.second[m_iDr0];
    const double volumeScaling_ij = con_i.second[m_iVolumeScaling];
    const double vol = vol_j * volumeScaling_ij;
    const double w = weightFunction(dr0);

    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
    }

    for (int d = 0; d < m_dim; d++) {
      for (int d2 = 0; d2 < m_dim; d2++) {
        K(d, d2) += w * dr0_ij[d] * dr0_ij[d2] * vol;
      }
    }
    nConnected++;
  }

  if (nConnected <= 3) {
    if (m_dim >= 2) {
      m_data(i, m_indexK[0]) = 0;
      m_data(i, m_indexK[1]) = 0;
      m_data(i, m_indexK[2]) = 0;
      m_data(i, m_indexK[3]) = 0;
    }
    if (m_dim == 3) {
      m_data(i, m_indexK[3]) = 0;
      m_data(i, m_indexK[4]) = 0;
      m_data(i, m_indexK[5]) = 0;
    }
  } else {
    K = inv(K);
  }

  if (m_dim >= 2) {
    m_data(i, m_indexK[0]) = K(0, 0);
    m_data(i, m_indexK[1]) = K(1, 1);
    m_data(i, m_indexK[2]) = K(0, 1);
    m_data(i, m_indexK[3]) = K(1, 0);
  }
  if (m_dim == 3) {
    m_data(i, m_indexK[3]) = K(2, 2);
    m_data(i, m_indexK[4]) = K(0, 2);
    m_data(i, m_indexK[5]) = K(1, 2);
  }
}
//------------------------------------------------------------------------------
double PD_NOPD::calculateStableMass(const int id_a, const int a, double dt) {
  // TMP: using the stable mass of PMB-force for testing
  dt *= 1.1;
  const arma::mat &R0 = m_particles.r0();
  double m[m_dim];
  double dr0[m_dim];
  double k[m_dim];
  //--------------------------------------------------------------------------
  double nu;
  double k_ = 0;
  double c;

  if (m_dim == 3) {
    nu = 1. / 4.;
    k_ = m_E / (3. * (1. - 2. * nu));
    c = 18. * k_ / (M_PI * pow(m_delta, 4));
  } else if (m_dim == 2) {
    nu = 1. / 3.;
    k_ = m_E / (2. * (1. - nu));
    c = 12 * k_ / (m_h * M_PI * pow(m_delta, 3));
  } else if (m_dim == 1) {
    nu = 1. / 4.;
    k_ = m_E;
    c = 2 * m_E / (m_h * m_h * pow(m_delta, 2));
  }

  double m_a = 0;
  const vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_a);

  for (auto &con : PDconnections) {
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_b = con.first;
    const int b = m_idToCol.at(id_b);

    double dr0Len = 0;

    for (int d = 0; d < m_dim; d++) {
      dr0[d] = R0(a, d) - R0(b, d);
      dr0Len += dr0[d] * dr0[d];
    }

    dr0Len = sqrt(dr0Len);
    const double vol_b = m_data(b, m_iVolume);
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double V = vol_b * volumeScaling;
    const double w = weightFunction(dr0Len);
    m_a += dr0Len * w * V;
  }

  c = 2. * k_ * pow(m_dim, 2.) / m_a;
  //--------------------------------------------------------------------------

  for (int d = 0; d < m_dim; d++) {
    m[d] = 0;
  }

  for (int i = 0; i < m_dim; i++) {
    for (int d = 0; d < m_dim; d++) {
      k[d] = 0;
    }
    for (auto &con : PDconnections) {
      if (con.second[m_iConnected] <= 0.5)
        continue;

      const int id_b = con.first;
      const int b = m_idToCol.at(id_b);

      double dr0Len = 0;

      for (int d = 0; d < m_dim; d++) {
        dr0[d] = R0(a, d) - R0(b, d);
        dr0Len += dr0[d] * dr0[d];
      }

      dr0Len = sqrt(dr0Len);
      const double vol_b = m_data(b, m_iVolume);
      const double volumeScaling = con.second[m_iVolumeScaling];
      const double V = vol_b * volumeScaling;

      const double dr0Len3 = pow(dr0Len, 3);
      const double w = weightFunction(dr0Len);
      const double coeff = c * w * V / dr0Len3;

      double sum = 0;

      for (int d = 0; d < m_dim; d++) {
        sum += fabs(dr0[d]);
      }
      sum *= fabs(dr0[i]) * coeff;

      k[i] += sum;
    }

    m[i] = k[i];
  }

  double stiffness = 0;

  for (int d = 0; d < m_dim; d++) {
    if (m[d] > stiffness) {
      stiffness = m[d];
    }
  }

  return 4. * 0.25 * pow(dt, 2) * stiffness;
  //    return 8.*0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
}
