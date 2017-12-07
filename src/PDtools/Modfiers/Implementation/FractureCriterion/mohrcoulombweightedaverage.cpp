#include "mohrcoulombweightedaverage.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
MohrCoulombWeightedAverage::MohrCoulombWeightedAverage(double mu, double C,
                                                       double T)
    : m_C(C), m_T(T) {
  m_phi = mu * M_PI / 180.;
  m_d = tan(m_phi);
  m_neededProperties = {pair<string, int>("stress", 1)};

  m_weights[0] = 1.;
  m_weights[1] = 50.;
  double tot = m_weights[0] + m_weights[1];
  m_weights[0] /= tot;
  m_weights[1] /= tot;

  m_hasStepOne = true;
}
//------------------------------------------------------------------------------
void MohrCoulombWeightedAverage::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_idToCol = &m_particles->getIdToCol_v();

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
}
//------------------------------------------------------------------------------
void MohrCoulombWeightedAverage::initialize() {}
//------------------------------------------------------------------------------
void MohrCoulombWeightedAverage::evaluateStepOne(const int id_i, const int i) {
  // First calculating the total stress on an node. The stress is averaged
  // all its bonds, and the stress state on a bond is the mean of the
  // stress at each material point.

  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  const mat &data = *m_data;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);
  double cos_theta = cos(M_PI / 2. + m_phi);
  double sin_theta = sin(M_PI / 2. + m_phi);
  double tan_theta = tan(0.5 * (M_PI + m_phi));
  tan_theta *= tan_theta;

  if (m_dim == 2) {
    const double sx_i = data(i, m_indexStress[0]);
    const double sy_i = data(i, m_indexStress[1]);
    const double sxy_i = data(i, m_indexStress[2]);
    const double tr1 = fabs(sx_i) + fabs(sy_i);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol)[id_j];

      if ((*m_data)(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      const double sx_j = data(j, m_indexStress[0]);
      const double sy_j = data(j, m_indexStress[1]);
      const double sxy_j = data(j, m_indexStress[2]);
      const double tr2 = fabs(sx_j) + fabs(sy_j);

      double sx, sy, sxy;
      if (tr1 > tr2) {
        sx = m_weights[0] * sx_j + m_weights[1] * sx_i;
        sy = m_weights[0] * sy_j + m_weights[1] * sy_i;
        sxy = m_weights[0] * sxy_j + m_weights[1] * sxy_i;
      } else {
        sx = m_weights[0] * sx_i + m_weights[1] * sx_j;
        sy = m_weights[0] * sy_i + m_weights[1] * sy_j;
        sxy = m_weights[0] * sxy_i + m_weights[1] * sxy_j;
      }

      const double first = 0.5 * (sx + sy);
      const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);
      const double p_1 = first - second;
      const double p_2 = first + second;

      const double shear = fabs(0.5 * (p_1 - p_2) * sin_theta);
      const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * cos_theta;

      if (shear >= fabs(m_C - m_d * normal) && normal < 0) {
        con.second[m_indexConnected] = 0;
      } else if (p_2 >= m_T && normal > 0) {
        con.second[m_indexConnected] = 0;
      }
    }
  }
}
//------------------------------------------------------------------------------
}
