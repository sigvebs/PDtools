#include "bondenergyfracture.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"

namespace PDtools {
//------------------------------------------------------------------------------
BondEnergyFracture::BondEnergyFracture(double delta, double G0,
                                       vector<Force *> &forces, double h)
    : m_forces(forces), m_delta(delta), m_G(G0), m_h(h) {
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void BondEnergyFracture::initialize() {
  if (m_dim == 3) {
    m_wc = 4. * m_G / (M_PI * pow(m_delta, 4));
  } else {
    m_wc = 9. * m_G / (4. * m_h * pow(m_delta, 3));
  }

  m_indexUnbreakable = m_particles->getParamId("unbreakable");
  m_indexStretch = m_particles->getPdParamId("stretch");
  m_indexDr0 = m_particles->getPdParamId("dr0");
  m_indexConnected = m_particles->getPdParamId("connected");
  m_indexForceScaling = m_particles->getPdParamId("forceScalingBond");
  m_indexMicromodulus = m_particles->getParamId("micromodulus");
  m_idToCol = &m_particles->idToCol();
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_data = &m_particles->data();

  m_state = false;
  m_broken = false;
  for (Force *force : m_forces) {
    force->updateState();
  }
}
//------------------------------------------------------------------------------
void BondEnergyFracture::evaluateStepOne(const int id_i, const int i) {
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;

  const double c_i = (*m_data)(i, m_indexMicromodulus);

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = (*m_idToCol)[id_j];

    if ((*m_data)(j, m_indexUnbreakable) >= 1)
      continue;
    if (con.second[m_indexConnected] <= 0.5)
      continue;

    const double c_j = (*m_data)(j, m_indexMicromodulus);
    const double g_ij = con.second[m_indexForceScaling];
    const double c = 0.5 * (c_i + c_j) * g_ij;
    const double s = con.second[m_indexStretch];
    const double dr0 = con.second[m_indexDr0];
    const double w = 0.5 * c * s * s * dr0;
    if (w > m_wc) {
      con.second[m_indexConnected] = 0;
      m_broken = true;
      (*m_data)(i, m_indexBrokenNow) = 1;
    }
  }
}
//------------------------------------------------------------------------------
void BondEnergyFracture::evaluateStepTwo() {
  if (m_broken) {
    m_state = true;
  } else {
    m_state = false;
  }
  m_broken = false;
}
//------------------------------------------------------------------------------
}
