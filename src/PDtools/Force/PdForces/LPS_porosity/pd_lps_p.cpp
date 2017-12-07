#include "pd_lps_p.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPS_POROSITY::PD_LPS_POROSITY(PD_Particles &particles, double phi_c,
                                 double n, bool planeStress, bool analyticalM)
    : Force(particles), m_planeStress(planeStress), m_phi_c(phi_c), m_n(n) {
  m_iMicromodulus = m_particles.registerParameter("micromodulus", 1);
  m_iTheta = m_particles.registerParameter("theta", 0);
  m_iThetaNew = m_particles.registerParameter("thetaNew", 0);
  m_iMass = m_particles.registerParameter("LPS_mass", 1);

  m_iVolume = m_particles.getParamId("volume");
  m_iPorosity = m_particles.getParamId("porosity");
  m_iRho = m_particles.getParamId("rho");
  m_iDr0 = m_particles.getPdParamId("dr0");
  m_iVolumeScaling = m_particles.getPdParamId("volumeScaling");
  m_iStretch = m_particles.registerPdParameter("stretch");
  m_iConnected = m_particles.getPdParamId("connected");
  m_indexBrokenNow = m_particles.registerParameter("brokenNow", 0);

  //    m_iK = m_particles.registerParameter("k_local");
  //    m_iE = m_particles.registerParameter("E_local");
  m_iA = m_particles.registerParameter("LPS_a");
  m_iB = m_particles.registerParameter("LPS_b");

  m_ghostParameters.push_back("volume");
  m_ghostParameters.push_back("theta");
  m_ghostParameters.push_back("LPS_mass");
  m_ghostParameters.push_back("LPS_a");
  m_ghostParameters.push_back("LPS_b");
  m_ghostParameters.push_back("micromodulus"); // For contact forces

  m_initialGhostParameters = {"volume", "theta", "LPS_mass", "LPS_a", "LPS_b"};
  m_hasUpdateState = true;
  m_analyticalM = analyticalM;

  //----------------------------------
  // TMP - standard stres calc from
  if (m_dim == 2) {
    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStress[1] = m_particles.registerParameter("s_yy");
    m_indexStress[2] = m_particles.registerParameter("s_xy");
  } else if (m_dim == 3) {
    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStress[1] = m_particles.registerParameter("s_yy");
    m_indexStress[2] = m_particles.registerParameter("s_xy");
    m_indexStress[3] = m_particles.registerParameter("s_zz");
    m_indexStress[4] = m_particles.registerParameter("s_xz");
    m_indexStress[5] = m_particles.registerParameter("s_yz");
  }
  //----------------------------------
  m_mass = m_data.colptr(m_iMass);
  m_theta = m_data.colptr(m_iTheta);
  m_volume = m_data.colptr(m_iVolume);
  m_porosity = m_data.colptr(m_iPorosity);
  m_alpha = m_data.colptr(m_iA);
  m_beta = m_data.colptr(m_iB);
  m_x = m_r.colptr(0);
  m_y = m_r.colptr(1);
  m_z = m_r.colptr(2);

  m_Fx = m_F.colptr(0);
  m_Fy = m_F.colptr(1);
  m_Fz = m_F.colptr(2);
  m_theta_new = m_data.colptr(m_iThetaNew);
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::calculateForces(const int id, const int i) {
  const double alpha_i = m_alpha[i];
  const double beta_i = m_beta[i];
  const double theta_i = m_theta[i];
  const double m_i = m_mass[i];

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);
  const int nConnections = PDconnections.size();
  double dr_ij[m_dim];

  double thetaNew = 0;
  int nConnected = 0;

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
      vector<double> &con_data = con.second;
      if (con_data[m_iConnected] <= 0.5)
        continue;

      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double m_j = m_mass[j];
      const double theta_j = m_theta[j];
      const double vol_j = m_volume[j];
      const double alpha_j = m_alpha[j];
      const double beta_j = m_beta[j];

      const double dr0 = con_data[m_iDr0];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      dr_ij[2] = m_z[j] - m_z[i];
      double dr2 =
          dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      const double dr = sqrt(dr2);
      const double ds = dr - dr0;
      double bond = (beta_i * theta_i / m_i + beta_j * theta_j / m_j) * dr0;
      bond += (alpha_i / m_i + alpha_j / m_j) * ds;
      bond *= w * vol_j * volumeScaling / dr;
      thetaNew += w * dr0 * ds * vol_j * volumeScaling;

      m_Fx[i] += dr_ij[0] * bond;
      m_Fy[i] += dr_ij[1] * bond;
      m_Fz[i] += dr_ij[2] * bond;

      nConnected++;

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
      const double alpha_j = m_alpha[j];
      const double beta_j = m_beta[j];

      const double dr0 = con_data[m_iDr0];
      const double volumeScaling = con_data[m_iVolumeScaling];
      const double w = weightFunction(dr0);

      dr_ij[0] = m_x[j] - m_x[i];
      dr_ij[1] = m_y[j] - m_y[i];
      double dr2 = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1];

      const double dr = sqrt(dr2);
      const double ds = dr - dr0;
      double bond = (beta_i * theta_i / m_i + beta_j * theta_j / m_j) * dr0;
      bond += (alpha_i / m_i + alpha_j / m_j) * ds;
      bond *= w * vol_j * volumeScaling / dr;
      thetaNew += w * dr0 * ds * vol_j * volumeScaling;

      m_Fx[i] += dr_ij[0] * bond;
      m_Fy[i] += dr_ij[1] * bond;

      nConnected++;
      //----------------------------------
      // TMP - standard stres calc from MD
      m_data(i, m_indexStress[0]) += 0.5 * dr_ij[0] * dr_ij[0] * bond;
      m_data(i, m_indexStress[1]) += 0.5 * dr_ij[1] * dr_ij[1] * bond;
      m_data(i, m_indexStress[2]) += 0.5 * dr_ij[0] * dr_ij[1] * bond;
    }
  }
  if (nConnected <= 3)
    m_theta_new[i] = 0;
  else
    m_theta_new[i] = m_t * thetaNew / m_i;

  m_continueState = false;
}
//------------------------------------------------------------------------------
double PD_LPS_POROSITY::calculatePotentialEnergyDensity(const int id_i,
                                                        const int i) {
  const double theta_i = this->computeDilation(id_i, i);
  double dr_ij[m_dim];
  const double a_i = m_data(i, m_iA);
  const double k_i = 1.; // TODO: TMP SOLUTION

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int nConnections = PDconnections.size();

  double W_i = 0;
  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double w = weightFunction(dr0);
    const double volumeScaling = con.second[m_iVolumeScaling];
    double dr2 = 0;

    for (int d = 0; d < m_dim; d++) {
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    const double ds = dr - dr0;
    const double extension_term =
        a_i * w * (pow(ds - theta_i * dr0 / m_dim, 2));

    W_i += (extension_term)*vol_j * volumeScaling;
  }
  W_i += k_i * (pow(theta_i, 2));

  return 0.5 * W_i;
}
//------------------------------------------------------------------------------
double PD_LPS_POROSITY::computeDilation(const int id_i, const int i) {
  const double m_i = m_data(i, m_iMass);
  double dr_ij[m_dim];

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int nConnections = PDconnections.size();

  double theta_i = 0;

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double volumeScaling = con.second[m_iVolumeScaling];
    double dr2 = 0;
    const double w = weightFunction(dr0);

    for (int d = 0; d < m_dim; d++) {
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    const double ds = dr - dr0;
    theta_i += w * dr0 * ds * vol_j * volumeScaling;
  }

  if (nConnections <= 3) {
    theta_i = 0;
  }
  return theta_i * m_t / m_i;
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::calculatePotentialEnergy(const int id_i, const int i,
                                               int indexPotential) {
  double vol_i = m_data(i, m_iVolume);
  m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i) * vol_i;
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::updateState(int id, int i) {
  (void)id;
  m_data(i, m_iTheta) = m_data(i, m_iThetaNew);
}
//------------------------------------------------------------------------------
double PD_LPS_POROSITY::calculateStableMass(const int id_a, const int a,
                                            double dt) {
  dt *= 1.1;

//  const double vol_a = m_data(a, m_iVolume);
  const double m_a = m_data(a, m_iMass);
  const arma::mat &R0 = m_particles.r0();
  const double alpha_a = m_data(a, m_iA);

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
      //            const double Va = vol_a*volumeScaling;
      const double Vb = vol_b * volumeScaling;
      const double w = weightFunction(dr0Len);
      const double alpha_b = m_data(b, m_iA);

      // Check this
      const double dr0Len2 = pow(dr0Len, 2);
      //            double C = m_dim*m_c*(Vb/pow(m_a,2) + Va/pow(m_b,2)) +
      //            m_alpha*(1./m_a + 1./m_b);
      double C = (alpha_a / m_a + alpha_b / m_b);
      //            double C = m_dim*m_c*(1./pow(m_a,2) + 1./pow(m_b,2)) +
      //            m_alpha*(1./m_a + 1./m_b);
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

  if (stiffness == 0)
    stiffness = 1.;

  return 4. * 0.25 * pow(dt, 2) * stiffness;
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::initialize(double E, double nu, double delta, int dim,
                                 double h, double lc) {
  Force::initialize(E, nu, delta, dim, h, lc);
  m_delta = delta;
  m_dim = dim;
  m_nu = nu;

  const int nParticles = m_particles.nParticles();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const double porosity = m_data(i, m_iPorosity);
    //        m_data(i, m_iVolume) *= (1-porosity);
    m_data(i, m_iRho) *= (1 - porosity);
  }

  for (int i = 0; i < nParticles; i++) {
    const double porosity = m_data(i, m_iPorosity);

    //        const double phi = 1 - porosity;
    const double E_loc = E * pow((1 - porosity / m_phi_c), m_n);
    const double nu_loc = nu * (1 - porosity / m_phi_c);
    //        const double nu_loc = nu;

    double mu = 0.5 * E_loc / (1 + nu_loc);
    double k = E_loc / (3. * (1. - 2. * nu_loc));

    double alpha, beta;

    if (dim == 3) {
      m_t = 3.;
      beta = (3. * k - 5. * mu);
      alpha = 15. * mu;
    } else if (dim == 2) {
      alpha = 8. * mu;

      if (m_planeStress) {
        m_t = 2. * (2. * nu_loc - 1.) / (nu_loc - 1.);
        k = k + mu / 9. * pow((nu_loc + 1.) / (2 * nu_loc - 1), 2);
        beta = m_t * k - 8. / 3. * mu * (2. - m_t / 3.);
      } else { // Plane strain
        m_t = 2.;
        beta = 2. * (k - 15. / 9. * mu);
        k = k + mu / 9.;
      }
    } else {
      cerr << "ERROR: dimension " << dim << " not supported." << endl;
      cerr << "use 2 or 3." << endl;
      exit(EXIT_FAILURE);
    }

    m_alpha[i] = alpha;
    m_beta[i] = beta;
  }

  calculateWeightedVolume();
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::calculateWeightedVolume() {
  const int nParticles = m_particles.nParticles();

// Calculating the one-body forces
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const int id_i = m_colToId.at(i);

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
      const double vol_j = m_data(j, m_iVolume);
      const double dr0 = con.second[m_iDr0];
      const double w = weightFunction(dr0);

      m += w * dr0 * dr0 * vol_j * volumeScaling;
    }

    // Remember to use the correct weightfunction
    if (m_analyticalM || nConnections < 3) {
      if (m_dim == 3) {
        //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
        m = M_PI * pow(m_delta, 5); //  For w = d/dr0
      } else {
        //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
        m = 2. * m_h * M_PI / 3. * pow(m_delta, 4); // For w = d/dr0
      }
    }

    m_data(i, m_iMass) = m;
    const double alpha = m_data(i, m_iA);
    m_data(i, m_iMicromodulus) = m_delta * alpha / m;
  }
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::updateWeightedVolume(int id_i, int i) {
  const vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int nConnections = PDconnections.size();

  double m = 0;
  int nActiveConnections = 0;
  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double w = weightFunction(dr0);

    m += w * dr0 * dr0 * vol_j * volumeScaling;
    nActiveConnections++;
  }

  if (m_analyticalM) {
    if (m_dim == 3) {
      //                m = 4.*M_PI/5.*pow(m_delta, 5); // For w = 1.
      m = M_PI * pow(m_delta, 5); //  For w = d/dr0
    } else {
      //                m = m_h*M_PI/2.*pow(m_delta, 4);// For w = 1.
      m = 2. * m_h * M_PI / 3. * pow(m_delta, 4); // For w = d/dr0
    }
  }

  if (nActiveConnections < 5) {
    m = m_data(i, m_iMass);
  }

  m_data(i, m_iMass) = m;

  // Setting the micromodulus
  const double alpha = m_data(i, m_iA);
  m_data(i, m_iMicromodulus) = m_delta * alpha / m;
}
//------------------------------------------------------------------------------
void PD_LPS_POROSITY::calculateStress(const int id_i, const int i,
                                      const int (&indexStress)[6]) {
  const double theta_i = m_data(i, m_iTheta);
  //    const double theta_i = this->computeDilation(id, i);
  const double m_i = m_data(i, m_iMass);
  const double a_i = m_data(i, m_iA);
  const double b_i = m_data(i, m_iB);

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);

  const int nConnections = PDconnections.size();
  double dr_ij[m_dim];
  double dr0_ij[m_dim];

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];
    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double m_j = m_data(j, m_iMass);
    const double theta_j = m_data(j, m_iTheta);
    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double w = weightFunction(dr0);
    const double a_j = m_data(j, m_iA);
    const double b_j = m_data(j, m_iB);

    double dr2 = 0;

    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    const double ds = dr - dr0;
    double bond = (b_i * theta_i / m_i + b_j * theta_j / m_j) * dr0;
    bond += (a_i / m_i + a_j / m_j) * ds;
    bond *= w * vol_j * volumeScaling / dr;

    m_data(i, indexStress[0]) += 0.5 * bond * dr_ij[X] * dr0_ij[X];
    m_data(i, indexStress[1]) += 0.5 * bond * dr_ij[Y] * dr0_ij[Y];
    m_data(i, indexStress[2]) += 0.5 * bond * dr_ij[X] * dr0_ij[Y];

    if (m_dim == 3) {
      m_data(i, indexStress[3]) += 0.5 * bond * dr_ij[Z] * dr0_ij[Z];
      m_data(i, indexStress[4]) += 0.5 * bond * dr_ij[X] * dr0_ij[Z];
      m_data(i, indexStress[5]) += 0.5 * bond * dr_ij[Y] * dr0_ij[Z];
    }
  }
}
//------------------------------------------------------------------------------
}
