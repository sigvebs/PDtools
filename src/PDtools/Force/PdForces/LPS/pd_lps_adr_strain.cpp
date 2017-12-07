#include "pd_lps_adr_strain.h"

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPS_ADR_STRAIN::PD_LPS_ADR_STRAIN(PD_Particles &particles, double c,
                                     double stretchCrit, double shearCrit,
                                     bool planeStress, bool analyticalM)
    : PD_LPS_CRIT_STRAIN(particles, c, stretchCrit, shearCrit, planeStress,
                         analyticalM) {}
//------------------------------------------------------------------------------
void PD_LPS_ADR_STRAIN::evaluateStatic(int id_i, int i) {
  if (m_data(i, m_iUnbreakable) >= 1)
    return;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  bool broken = false;
  const double theta_i = m_data(i, m_iTheta);

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = m_idToCol[id_j];

    if (m_data(j, m_iUnbreakable) >= 1)
      continue;

    if (con.second[m_iConnected] <= 0.5)
      continue;

    const double s = con.second[m_iStretch];
    const double theta_j = m_data(j, m_iTheta);
    const double s_d = std::max(s - theta_i / 3., s - theta_j / 3.);

    if (s > m_stretchCrit) {
      m_data(i, m_indexBrokenNow) = 1;
      con.second[m_iConnected] = 0;
      m_continueState = true;
      broken = true;
    } else if (s_d > m_shearCrit) {
      m_data(i, m_indexBrokenNow) = 1;
      con.second[m_iConnected] = 0;
      m_continueState = true;
      broken = true;
    }
  }

  if (broken) {
    updateWeightedVolume(id_i, i);
    //        m_data(i, m_indexBrokenNow) = 0;
  }
}
//------------------------------------------------------------------------------
void PD_LPS_ADR_STRAIN::evaluateStepTwo(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
}
