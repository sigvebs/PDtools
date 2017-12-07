#include "force.h"

#include "Particles/pd_particles.h"
#include <stdio.h>

namespace PDtools {
//------------------------------------------------------------------------------
vector<pair<string, int>> Force::getNeededProperties() const {
  return neededProperties;
}
//------------------------------------------------------------------------------
int Force::getCalulateStress() const { return m_calulateStress; }
//------------------------------------------------------------------------------
void Force::evaluateStepOne() {}
//------------------------------------------------------------------------------
bool Force::getHasStepOneModifier() const { return m_hasStepOneModifier; }
//------------------------------------------------------------------------------
void Force::evaluateStepOne(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
void Force::evaluateStepTwo(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
void Force::evaluateStatic(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
bool Force::getHasStepTwoModifier() const { return m_hasStepTwoModifier; }
//------------------------------------------------------------------------------
bool Force::getContinueState() const { return m_continueState; }
//------------------------------------------------------------------------------
bool Force::getHasStaticModifier() const { return m_hasStaticModifier; }
//------------------------------------------------------------------------------
bool Force::getHasUpdateState() const { return m_hasUpdateState; }
//------------------------------------------------------------------------------
Force::Force(PD_Particles &particles, string _type)
    : m_particles(particles), m_r(m_particles.r()), m_v(m_particles.v()),
      m_r0(m_particles.r0()), m_F(m_particles.F()), m_data(m_particles.data()),
      m_idToCol_v(m_particles.getIdToCol_v()),
      m_colToId(m_particles.colToId()),
      m_idToElement(m_particles.getIdToElement()),
      m_triElements(m_particles.getTriElements()),
      m_quadElements(m_particles.getQuadElements()), m_dim(m_particles.dim()),
      name(_type) {}
//------------------------------------------------------------------------------
Force::~Force() {}
//------------------------------------------------------------------------------
void Force::initialize(double E, double nu, double delta, int dim, double h,
                       double lc) {
  m_E = E;
  m_nu = nu;
  m_h = h;
  m_delta = delta;
  m_dim = dim;
  m_lc = lc;
}
//------------------------------------------------------------------------------
double Force::calculatePotentialEnergyDensity(const int id_i, const int i) {
  (void)id_i;
  (void)i;
  return 0;
}
//------------------------------------------------------------------------------
void Force::calculatePotentialEnergy(const int id_i, const int i,
                                     int indexPotential) {
  (void)id_i;
  (void)i;
  (void)indexPotential;
  //    cerr << "ERROR: potential energy not implemented for this force" <<
  //    endl;
}
//------------------------------------------------------------------------------
double Force::calculateBondEnergy(const pair<int, int> &idCol,
                                  pair<int, vector<double>> &con) {
  (void)idCol;
  (void)con;

  return 0;
}
//------------------------------------------------------------------------------
void Force::calculateStress(const int id_i, const int i,
                            const int (&indexStress)[6]) {
  (void)id_i;
  (void)i;
  (void)indexStress;
  //    cerr << "ERROR: stress not implemented for this force" << endl;
}
//------------------------------------------------------------------------------
void Force::updateState() { m_counter++; }
//------------------------------------------------------------------------------
void Force::updateState(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
double Force::calculateStableMass(const int id_i, const int i, double dt) {
  (void)id_i;
  (void)i;
  (void)dt;

  return 0.;
}
//------------------------------------------------------------------------------
void Force::numericalInitialization(bool ni) { m_numericalInitialization = ni; }
//------------------------------------------------------------------------------
void Force::setDim(int dim) { m_dim = dim; }
//------------------------------------------------------------------------------
int Force::initializeSurfaceCorrection() {
  if (!m_hasSurfaceCorrection)
    return 0;

  // Temporary registering the correction vectors
  for (int d = 0; d < m_dim; d++) {
    m_g(d) = m_particles.registerParameter(forceScalingStringIds[d]);
  }

  return 1;
}
//------------------------------------------------------------------------------
void Force::applySurfaceCorrectionStep1(double strain) {
  applyStrainCorrection(strain);
  //    applyShearCorrection(strain);
}
//------------------------------------------------------------------------------
void Force::applyStrainCorrection(double strain) {
  const ivec &colToId = m_particles.colToId();
  const size_t nParticles = m_particles.nParticles();
  const size_t nParticlesAndGhosts =
      m_particles.nParticles() + m_particles.nGhostParticles();

  // Stretching all particle in the x-direction
  arma::vec3 scaleFactor;
  scaleFactor(0) = strain;
  scaleFactor(1) = 0;
  scaleFactor(2) = 0;

  double W_bulk = 0; // Analytical bulk energy
  switch (m_dim) {
  case 3:
    W_bulk = 0.6 * m_E * pow(strain, 2);
    break;
  case 2:
    W_bulk = 9. / 16. * m_E * pow(strain, 2);
    break;
  case 1:
    W_bulk = 0.5 * m_E * pow(strain, 2);
    break;
  }

  for (int a = 0; a < m_dim; a++) {
    if (a == 1)
      scaleFactor.swap_rows(0, 1);
    else if (a == 2)
      scaleFactor.swap_rows(1, 2);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Loading the geometry
    for (size_t i = 0; i < nParticlesAndGhosts; i++) {
      for (int d = 0; d < m_dim; d++) {
        m_r(i, d) = (1 + scaleFactor(d)) * m_r(i, d);
      }
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Calculating the elastic energy density
    for (int i = 0; i < nParticles; i++) {
      const int id_i = colToId(i);
      const double W = this->calculatePotentialEnergyDensity(id_i, i);
      m_data(i, m_g(a)) = W_bulk / W;
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Resetting the positions
    for (size_t i = 0; i < nParticlesAndGhosts; i++) {
      for (int d = 0; d < m_dim; d++) {
        m_r(i, d) = m_r(i, d) / (1 + scaleFactor(d));
      }
    }
  }
}
//------------------------------------------------------------------------------
void Force::applySurfaceCorrectionStep2() {
  const ivec &colToId = m_particles.colToId();
  const int iDr0 = m_particles.getPdParamId("dr0");
  const int iForceScaling = m_particles.getPdParamId("forceScalingBond");
  const int nParticles = m_particles.nParticles();

// Calculating the scaling
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const int id_i = colToId(i);
    vector<pair<int, vector<double>>> &PDconnections =
        m_particles.pdConnections(id_i);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      const double dr0Len = con.second[iDr0];
      const arma::vec3 &n = (m_r.row(i).t() - m_r.row(j).t()) / dr0Len;

      arma::vec3 g_mean;
      double G = 0;
      for (int d = 0; d < m_dim; d++) {
        const double g_i = m_data(i, m_g(d));
        const double g_j = m_data(j, m_g(d));
        g_mean(d) = 0.5 * (g_i + g_j);
        G += pow(n(d) / g_mean(d), 2);
      }
      G = pow(G, -0.5);
      con.second[iForceScaling] *= G;
    }
  }
}
//------------------------------------------------------------------------------
void Force::applyShearCorrection(double shear) {
  const ivec &colToId = m_particles.colToId();
  double m_mu = 0.5 * m_E / (1 + m_nu); // Tmp
  arma::vec3 strainFactor;

  // Performing a simple shear of all particle in the x, y and z-direction
  arma::ivec3 axis;
  strainFactor(0) = shear;
  strainFactor(1) = 0;
  strainFactor(2) = 0;
  axis(0) = 1;
  axis(1) = 0;
  axis(2) = 0;

  const double W_s = 0.5 * m_mu * shear * shear; // Analytical solution

  for (int a = 0; a < m_dim; a++) {
    if (a == 1) {
      strainFactor.swap_rows(1, 2);
      strainFactor.swap_rows(0, 1);
      axis(0) = 2;
      axis(1) = 0;
      axis(2) = 1;
    } else if (a == 2) {
      strainFactor.swap_rows(2, 0);
      strainFactor.swap_rows(1, 2);
      axis(0) = 2;
      axis(1) = 0;
      axis(2) = 0;
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Loading the geometry
    for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
      for (int d = 0; d < m_dim; d++) {
        const double shearFactor = strainFactor(d) * m_r(i, axis(d));
        m_r(i, d) = m_r(i, d) + shearFactor;
      }
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Calculating the elastic energy density
    for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
      const int id_i = colToId(i);
      double W = this->calculatePotentialEnergyDensity(id_i, i);
      m_data(i, m_g(a)) = W_s / W;
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Resetting the positions
    for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
      for (int d = 0; d < m_dim; d++) {
        m_r(i, d) = m_r0(i, d); // Tmp
      }
    }
  }
}
//------------------------------------------------------------------------------
vector<string> Force::initalGhostDependencies() {
  return m_initialGhostParameters;
}
//------------------------------------------------------------------------------
vector<string> Force::ghostDependencies() { return m_ghostParameters; }
//------------------------------------------------------------------------------
vector<string> Force::getSurfaceCorrectionGhostParameters() {
  vector<string> param;
  for (int d = 0; d < m_dim; d++) {
    param.push_back(forceScalingStringIds[d]);
  }
  return param;
}
//------------------------------------------------------------------------------
}
