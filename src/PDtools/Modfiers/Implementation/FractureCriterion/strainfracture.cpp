#include "strainfracture.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
StrainFracture::StrainFracture(double Eeq, double Evol, int dim)
    : m_Eeq(Eeq), m_Evol(Evol), m_dim(dim) {
  m_neededProperties = {pair<string, int>("strain", 1)};
  m_hasStepOne = true;
}
//------------------------------------------------------------------------------
void StrainFracture::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexCompute = m_particles->registerPdParameter("compute");
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_particles->registerParameter("damage");
  m_idToCol = &m_particles->idToCol();

  switch (m_dim) {
  case 1:
    m_ghostParameters = {"e_xx"};
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    break;
  case 2:
    m_ghostParameters = {"e_xx", "e_yy", "e_xy"};
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    m_indexStrain[1] = m_particles->registerParameter("e_yy");
    m_indexStrain[2] = m_particles->registerParameter("e_xy");
    break;
  case 3:
    m_ghostParameters = {"e_xx", "e_yy", "e_zz", "e_xy", "e_xz", "e_yz"};
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    m_indexStrain[1] = m_particles->registerParameter("e_yy");
    m_indexStrain[2] = m_particles->registerParameter("e_xy");
    m_indexStrain[3] = m_particles->registerParameter("e_zz");
    m_indexStrain[4] = m_particles->registerParameter("e_xz");
    m_indexStrain[5] = m_particles->registerParameter("e_yz");
    break;
  }
}
//------------------------------------------------------------------------------
void StrainFracture::evaluateStepOne(const int id_i, const int i) {
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

      const double ex =
          0.5 * (data(i, m_indexStrain[0]) + data(j, m_indexStrain[0]));
      const double ey =
          0.5 * (data(i, m_indexStrain[1]) + data(j, m_indexStrain[1]));
      const double exy =
          0.5 * (data(i, m_indexStrain[2]) + data(j, m_indexStrain[2]));

      const double Eeq = 2. / 9. * sqrt(pow(ex - ey, 2) + ex * ex + ey * ey +
                                        4. / 3. * exy * exy);
      const double Evol = ex + ey + ex * ey;

      if (Eeq >= m_Eeq) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      } else if (Evol >= m_Evol) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }
    }
  } else if (m_dim == 3) {
  }
}
//------------------------------------------------------------------------------
}
