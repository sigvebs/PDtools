#include "pd_bondforcegaussian.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_bondforceGaussian::PD_bondforceGaussian(PD_Particles &particles,
                                           std::string weightType)
    : Force(particles), m_weightType(weightType) {
  m_indexMicromodulus = m_particles.registerParameter("micromodulus");
  m_indexVolume = m_particles.getParamId("volume");
  m_indexDr0 = m_particles.getPdParamId("dr0");
  m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
  m_indexForceScaling = m_particles.registerPdParameter("forceScalingBond", 1);
  m_indexStretch = m_particles.registerPdParameter("stretch");
  m_indexConnected = m_particles.getPdParamId("connected");
  m_indexWeightFunction = m_particles.registerPdParameter("weightFunction", 1);

  m_hasSurfaceCorrection = true;
  m_ghostParameters = {"volume", "micromodulus"};
  m_initialGhostParameters = {"volume", "micromodulus"};
}
//------------------------------------------------------------------------------
PD_bondforceGaussian::~PD_bondforceGaussian() {}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculateForces(const int id, const int i) {
  // PD_bond
  const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
  const double vol_i = m_data(i, m_indexVolume);
  const int nParticles = m_particles.nParticles();
#endif
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);

  double dr_ij[m_dim];

  for (auto &con : PDconnections) {
    if (con.second[m_indexConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

#if USE_N3L
    if (j < i)
      continue;
#endif
    const double c_j = m_data(j, m_indexMicromodulus);
    const double vol_j = m_data(j, m_indexVolume);
    const double dr0 = con.second[m_indexDr0];
    const double volumeScaling = con.second[m_indexVolumeScaling];
    const double g_ij = con.second[m_indexForceScaling];
    const double w_ij = con.second[m_indexWeightFunction];
    const double c_ij = 0.5 * (c_i + c_j) * w_ij * g_ij;

    double dr2 = 0;
    for (int d = 0; d < m_dim; d++) {
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    double ds = dr - dr0;

    // To avoid roundoff errors
    if (fabs(ds) < THRESHOLD)
      ds = 0.0;

    const double s = ds / dr0;
    const double fbond_ij = c_ij * s * vol_j * volumeScaling / dr;

    for (int d = 0; d < m_dim; d++) {
      m_F(i, d) += dr_ij[d] * fbond_ij;
    }

    con.second[m_indexStretch] = s;
#if USE_N3L
    if (j > i && j < nParticles) {
      const double fbond_ji = -c_ij * s * vol_i * volumeScaling / dr;
      for (int d = 0; d < m_dim; d++) {
        m_F(j, d) += dr_ij[d] * fbond_ji;
      }
    }
#endif
  }
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculatePotentialEnergyDensity(const int id_i,
                                                             const int i) {
  // PD_bond
  const double c_i = m_data(i, m_indexMicromodulus);
  double dr_ij[m_dim];
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);

  double energy = 0;
  for (auto &con : PDconnections) {
    if (con.second[m_indexConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double vol_j = m_data(j, m_indexVolume);
    const double dr0 = con.second[m_indexDr0];
    const double volumeScaling = con.second[m_indexVolumeScaling];
    const double c_j = m_data(j, m_indexMicromodulus);
    const double g_ij = con.second[m_indexForceScaling];
    const double w_ij = con.second[m_indexWeightFunction];
    const double c_ij = 0.5 * (c_i + c_j) * w_ij * g_ij;

    double dr2 = 0;
    for (int d = 0; d < m_dim; d++) {
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    double s = (dr - dr0) / dr0;
    energy += c_ij * (s * s) * dr0 * vol_j * volumeScaling;
  }

  return 0.25 * energy;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculatePotentialEnergy(const int id_i, const int i,
                                                    int indexPotential) {
  // PD_bond
  double vol_i = m_data(i, m_indexVolume);
  m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i) * vol_i;
}
//------------------------------------------------------------------------------
double
PD_bondforceGaussian::calculateBondEnergy(const int id_i, const int i,
                                          pair<int, vector<double>> &con) {
  // PD_bond
  (void)id_i;
  const double c_i = m_data(i, m_indexMicromodulus);
  const int id_j = con.first;
  const int j = m_idToCol_v[id_j];

  const double vol_j = m_data(j, m_indexVolume);
  const double dr0Len = con.second[m_indexDr0];
  const double volumeScaling = con.second[m_indexVolumeScaling];
  const double g_ij = con.second[m_indexForceScaling];
  const double w_ij = con.second[m_indexWeightFunction];
  const double c_j = m_data(j, m_indexMicromodulus);
  const double c = 0.5 * (c_i + c_j) * w_ij * g_ij;

  double dr_ij[m_dim];
  double dr2 = 0;

  for (int d = 0; d < m_dim; d++) {
    dr_ij[d] = m_r(j, d) - m_r(i, d);
    dr2 += dr_ij[d] * dr_ij[d];
  }

  const double drLen = sqrt(dr2);
  double ds = drLen - dr0Len;

  // To avoid roundoff errors
  if (fabs(ds) < THRESHOLD)
    ds = 0.0;

  const double energy = 0.5 * c * (ds * ds) / dr0Len * vol_j * volumeScaling;

  return energy;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculateStress(const int id_i, const int i,
                                           const int (&indexStress)[6]) {
  // PD_bond
  const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
  const double vol_i = m_data(i, m_indexVolume);
  const int nParticles = m_particles.nParticles();
#endif

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  double dr_ij[m_dim];
  double f[m_dim];

  for (auto &con : PDconnections) {
    if (con.second[m_indexConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double c_j = m_data(j, m_indexMicromodulus);
    const double vol_j = m_data(j, m_indexVolume);
    const double dr0 = con.second[m_indexDr0];
    const double volumeScaling = con.second[m_indexVolumeScaling];
    const double g_ij = con.second[m_indexForceScaling];
    const double w_ij = con.second[m_indexWeightFunction];
    const double c_ij = 0.5 * (c_i + c_j) * w_ij * g_ij;
    double dr2 = 0;

    for (int d = 0; d < m_dim; d++) {
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    double ds = dr - dr0;

    // To avoid roundoff errors
    if (fabs(ds) < THRESHOLD)
      ds = 0.0;

    double s = ds / dr0;
    //        double stretch = con.second[m_indexStretch];
    const double fbond_ij = c_ij * s * vol_j * volumeScaling / dr;

    for (int d = 0; d < m_dim; d++) {
      f[d] = dr_ij[d] * fbond_ij;
    }
    m_data(i, indexStress[0]) += 0.5 * f[X] * dr_ij[X];
    m_data(i, indexStress[1]) += 0.5 * f[Y] * dr_ij[Y];
    m_data(i, indexStress[2]) += 0.5 * f[X] * dr_ij[Y];

    if (m_dim == 3) {
      m_data(i, indexStress[3]) += 0.5 * f[Z] * dr_ij[Z];
      m_data(i, indexStress[4]) += 0.5 * f[X] * dr_ij[Z];
      m_data(i, indexStress[5]) += 0.5 * f[Z] * dr_ij[Z];
    }

#if USE_N3L
    if (j > i && j < nParticles) {
      const double bond_ji = c_ij * s * vol_i * volumeScaling / dr;
      m_data(j, indexStress[0]) += 0.5 * bond_ji * dr_ij[X] * dr_ij[X];
      m_data(j, indexStress[1]) += 0.5 * bond_ji * dr_ij[Y] * dr_ij[Y];
      m_data(j, indexStress[2]) += 0.5 * bond_ji * dr_ij[X] * dr_ij[Y];

      if (m_dim == 3) {
        m_data(j, indexStress[3]) += 0.5 * bond_ji * dr_ij[Z] * dr_ij[Z];
        m_data(j, indexStress[4]) += 0.5 * bond_ji * dr_ij[X] * dr_ij[Z];
        m_data(j, indexStress[5]) += 0.5 * bond_ji * dr_ij[Y] * dr_ij[Z];
      }
    }
#endif
  }
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculateStableMass(const int id_a, const int a,
                                                 double dt) {
  // PD_bond
  dt *= 1.1;
  const double c_a = m_data(a, m_indexMicromodulus);

  double m[m_dim];
  double dR0[m_dim];

  const arma::mat &matR0 = m_particles.r0();
  for (int d = 0; d < m_dim; d++) {
    m[d] = 0;
  }

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_a);

  double k[m_dim];

  for (int i = 0; i < m_dim; i++) {
    for (int d = 0; d < m_dim; d++) {
      k[d] = 0;
    }

    for (auto &con : PDconnections) {
      if (con.second[m_indexConnected] <= 0.5)
        continue;

      const int id_b = con.first;
      const int b = m_idToCol_v[id_b];

      double sum = 0;
      for (int d = 0; d < m_dim; d++) {
        dR0[d] = matR0(a, d) - matR0(b, d);
        sum += fabs(dR0[d]);
      }

      const double c_b = m_data(b, m_indexMicromodulus);
      const double vol_b = m_data(b, m_indexVolume);
      const double dr0 = con.second[m_indexDr0];
      const double volumeScaling = con.second[m_indexVolumeScaling];
      const double g_ab = con.second[m_indexForceScaling];
      const double w_ij = con.second[m_indexWeightFunction];
      const double c_ab = 0.5 * (c_a + c_b) * w_ij * g_ab;
      const double dr0_3 = pow(dr0, 3);

      const double coeff = c_ab * volumeScaling * vol_b / dr0_3;
      k[i] += coeff * fabs(dR0[i]) * sum;
    }

    m[i] = k[i];
  }

  double stiffness = 0;

  for (int d = 0; d < m_dim; d++) {
    if (m[d] > stiffness) {
      stiffness = m[d];
    }
  }

  return 2 * 0.25 * pow(dt, 2) * stiffness;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initialize(double E, double nu, double delta,
                                      int dim, double h, double lc) {
  Force::initialize(E, nu, delta, dim, h, lc);

  if (m_weightType == "constant") {
    initializeConstant();
  } else if (m_weightType == "linear") {
    initializeLinear();
  } else if (m_weightType == "gaussian") {
    initializeGaussian();
  } else if (m_weightType == "sigmoid") {
    initializeSigmoid();
  } else {
    std::cerr << "Weightfunction " << m_weightType << endl
              << "for bondforce" << std::endl;
    exit(EXIT_FAILURE);
  }
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeConstant() {
  double c, k;

  if (m_dim == 3) {
    m_nu = 1. / 4.;
    k = m_E / (3. * (1. - 2. * m_nu));
    c = 18. * k / (M_PI * pow(m_delta, 4));
  } else if (m_dim == 2) {
    m_nu = 1. / 3.;
    k = m_E / (2. * (1. - m_nu));
    c = 12 * k / (m_h * M_PI * pow(m_delta, 3));
  } else if (m_dim == 1) {
    m_nu = 1. / 4.;
    k = m_E;
    c = 2 * m_E / (m_h * m_h * pow(m_delta, 2));
  } else {
    cerr << "ERROR: dimension " << m_dim << " not supported" << endl;
    cerr << "use 1, 2 or 3." << endl;
    exit(EXIT_FAILURE);
  }

  m_particles.setParameter("micromodulus", c);
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeLinear() {
  double c = -1;
  double k = -1;

  if (m_dim == 3) {
    m_nu = 1. / 4.;
    k = m_E / (3. * (1. - 2. * m_nu));
    c = 90. * k / (M_PI * pow(m_delta, 4));
  } else if (m_dim == 2) {
    m_nu = 1. / 3.;
    k = m_E / (2. * (1. - m_nu));
    c = 48. * k / (m_h * M_PI * pow(m_delta, 3));
    c = 48. * k / (m_h * M_PI * pow(m_delta, 3) *
                   (4 - 3 * m_delta / (m_delta + 0.5 * m_lc)));
  } else if (m_dim == 1) {
    cerr << "dim == 1 is not implemented for linear weightfunction" << endl;
    exit(EXIT_FAILURE);
  }

  m_particles.setParameter("micromodulus", c);

  const ivec &colToId = m_particles.colToId();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
    const int pId = colToId(i);

    vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(pId);

    for (auto &con : PDconnections) {
      const double dr0 = con.second[m_indexDr0];
      //            con.second[m_indexWeightFunction] = 1 - dr0/m_delta;
      con.second[m_indexWeightFunction] = 1 - dr0 / (m_delta + 0.5 * m_lc);
    }
  }
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeGaussian() {
  double c = -1;
  double k = -1;
  m_l = 0.5 * m_delta;

  if (m_dim == 3) {
    m_nu = 1. / 4.;
    k = m_E / (3. * (1. - 2. * m_nu));
    c = 9. * k / (4 * M_PI * pow(m_delta, 4) * (3. - 8. * exp(-1)));
  } else if (m_dim == 2) {
    m_nu = 1. / 3.;
    k = m_E / (2. * (1. - m_nu));
    c = 4. * k / (m_h * M_PI * pow(m_delta, 3) * (2. - 5. * exp(-1)));
  } else if (m_dim == 1) {
    m_nu = 1. / 4.;
    k = m_E;
    c = 2 * k / (m_h * m_h * pow(m_delta, 2)) * (1. - 2. * exp(-1));
    c = 2 * k / (m_h * m_h * pow(m_delta, 2));
  }

  m_particles.setParameter("micromodulus", c);
  const ivec &colToId = m_particles.colToId();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
    const int pId = colToId(i);

    vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(pId);
    for (auto &con : PDconnections) {
      const double dr0Len = con.second[m_indexDr0];
      const double e_drl = exp(-dr0Len / m_l);

      con.second[m_indexWeightFunction] = e_drl;
    }
  }
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeSigmoid() {
  double c, k;
  double alpha = 0.5 * m_delta;
  double beta = m_delta / 10.;

  if (m_dim == 3) {
    m_nu = 1. / 4.;
    k = m_E / (3. * (1. - 2. * m_nu));
    c = 18. * k / (M_PI * pow(m_delta, 4));
  } else if (m_dim == 2) {
    m_nu = 1. / 3.;
    k = m_E / (2. * (1. - m_nu));
    c = 12 * k / (m_h * M_PI * pow(m_delta, 3));
  } else if (m_dim == 1) {
    m_nu = 1. / 4.;
    k = m_E;
    c = 2 * m_E / (m_h * m_h * pow(m_delta, 2));
  } else {
    cerr << "ERROR: dimension " << m_dim << " not supported" << endl;
    cerr << "use 1, 2 or 3." << endl;
    exit(EXIT_FAILURE);
  }

  m_particles.setParameter("micromodulus", c);
  const ivec &colToId = m_particles.colToId();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
    const int pId = colToId(i);

    vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(pId);
    for (auto &con : PDconnections) {
      const double dr0 = con.second[m_indexDr0];
      const double e_drl = 1. / (1. + exp((dr0 - alpha) / beta));

      con.second[m_indexWeightFunction] = e_drl;
    }
  }
}
//------------------------------------------------------------------------------
}
