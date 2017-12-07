#include "mohrcoulombmax.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
MohrCoulombMax::MohrCoulombMax(double mu, double C, double T) : m_C(C), m_T(T) {
  m_phi = mu * M_PI / 180.;
  m_d = tan(m_phi);
  m_neededProperties = {pair<string, int>("stress", 1)};

  m_weight1 = 0.95;
  m_weight2 = 1. - m_weight1;

  m_hasStepOne = true;
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void MohrCoulombMax::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexCompute = m_particles->registerPdParameter("compute");
  //    m_indexStressCenter = m_particles->registerPdParameter("stressIndex");
  m_indexBroken = m_particles->registerParameter("broken", 0);
  m_idToCol = &m_particles->getIdToCol_v();

  m_particles->registerParameter("damage");

  switch (m_dim) {
  case 1:
    m_ghostParameters = {"s_xx"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    break;
  case 2:
    m_ghostParameters = {"s_xx", "s_yy", "s_xy"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    break;
  case 3:
    m_ghostParameters = {"s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStress[3] = m_particles->registerParameter("s_zz");
    m_indexStress[4] = m_particles->registerParameter("s_xz");
    m_indexStress[5] = m_particles->registerParameter("s_yz");
    break;
  }
  m_ghostParameters.push_back("broken");
}
//------------------------------------------------------------------------------
void MohrCoulombMax::initialize() { }
//------------------------------------------------------------------------------
void MohrCoulombMax::evaluateStepOne(const int id_i, const int i) {
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  const mat &data = *m_data;

#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
  arma::vec eigval(m_dim);
#endif
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);
  double cos_theta = cos(M_PI / 2. + m_phi);
  double sin_theta = sin(M_PI / 2. + m_phi);

  // My stress
  double sx, sy, sxy;
  const int broken_i = data(i, m_indexBroken);

  double w_i, w_j;

  if (m_dim == 2) {
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol)[id_j];

      if ((*m_data)(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      const int broken_j = data(j, m_indexBroken);

      if (broken_i == 0 && broken_j == 0) {
        continue;
      }

      // Both are broken in the same type of fracture
      if (broken_i == broken_j) {
        con.second[m_indexConnected] = 0;
        continue;
      }

      // Adjusting the weight
      if (broken_i > broken_j) {
        w_i = m_weight1;
        w_j = m_weight2;
      } else {
        w_i = m_weight2;
        w_j = m_weight1;
      }

      // One is broken. Combine the two stresses and weight acordingly
      sx = w_i * data(i, m_indexStress[0]) + w_j * data(j, m_indexStress[0]);
      sy = w_i * data(i, m_indexStress[1]) + w_j * data(j, m_indexStress[1]);
      sxy = w_i * data(i, m_indexStress[2]) + w_j * data(j, m_indexStress[2]);

      const double first = 0.5 * (sx + sy);
      const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

      const double s1 = first + second;
      const double s2 = first - second;
      const double p_1 = std::min(s1, s2);
      const double p_2 = std::max(s1, s2);

      const double shear = std::fabs(0.5 * (p_1 - p_2) * sin_theta);
      const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * cos_theta;

      if (shear >= std::fabs(m_C - m_d * normal) && normal < 0) {
        con.second[m_indexConnected] = 0;
      } else if (p_2 >= m_T && normal > 0) {
        con.second[m_indexConnected] = 0;
      }
    }
  }
}
//------------------------------------------------------------------------------
void MohrCoulombMax::evaluateStepTwo(const int id_i, const int i) {
  (void)id_i;
  // First calculating the total stress on an material point. The stress is
  // averaged
  // all its bonds, and the stress state on a bond is the mean of the
  // stress at each material point.

  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  mat &data = *m_data;
  data(i, m_indexBroken) = 0;

  double cos_theta = cos(M_PI / 2. + m_phi);
  double sin_theta = sin(M_PI / 2. + m_phi);

  if (m_dim == 2) {
    double sx, sy, sxy;
    sx = data(i, m_indexStress[0]);
    sy = data(i, m_indexStress[1]);
    sxy = data(i, m_indexStress[2]);

    const double first = 0.5 * (sx + sy);
    const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

    const double p_2 = first + second;
    const double p_1 = first - second;

    const double shear = fabs(0.5 * (p_1 - p_2) * sin_theta);
    const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * cos_theta;

    if (shear >= fabs(m_C - m_d * normal) && normal < 0) {
      data(i, m_indexBroken) = 2;
    } else if (p_2 >= m_T) {
      data(i, m_indexBroken) = 1;
    }
  }
}
//------------------------------------------------------------------------------
void MohrCoulombMax::evaluateStepTwo() {
  if (m_broken) {
    m_state = true;
  } else {
    m_state = false;
  }

  m_broken = false;
}
//------------------------------------------------------------------------------
}
