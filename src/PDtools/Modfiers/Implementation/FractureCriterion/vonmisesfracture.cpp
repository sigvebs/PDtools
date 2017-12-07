#include "vonmisesfracture.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
VonMisesFracture::VonMisesFracture(double sigma_y)
    : m_sigma_y(sigma_y), m_sigma_y2(sigma_y * sigma_y) {
  m_hasStepOne = true;
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void VonMisesFracture::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexCompute = m_particles->registerPdParameter("compute");
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_particles->registerParameter("damage");
  m_idToCol = &m_particles->idToCol();

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
void VonMisesFracture::evaluateStepOne(const int id_i, const int i) {
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  mat &data = *m_data;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);

  if (m_dim == 2) {
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol)[id_j];

      if ((*m_data)(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      const double sx =
          0.5 * (data(i, m_indexStress[0]) + data(j, m_indexStress[0]));
      const double sy =
          0.5 * (data(i, m_indexStress[1]) + data(j, m_indexStress[1]));
      const double sxy =
          0.5 * (data(i, m_indexStress[2]) + data(j, m_indexStress[2]));

      const double s = sx * sx - sx * sy + sy * sy + 3 * sxy * sxy;

      if (s >= m_sigma_y2) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }
    }
  }
}
//------------------------------------------------------------------------------
}
