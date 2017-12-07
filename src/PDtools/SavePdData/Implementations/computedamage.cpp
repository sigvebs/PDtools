#include "computedamage.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
ComputeDamage::ComputeDamage(PD_Particles &particles)
    : ComputeProperty(particles) {
  m_indexDamage = m_particles.registerParameter("damage");
  m_indexMaxPdConnections = m_particles.registerParameter("maxPdConnections");
  m_data = &m_particles.data();
  m_indexConnected = m_particles.getPdParamId("connected");
}
//------------------------------------------------------------------------------
ComputeDamage::~ComputeDamage() {}
//------------------------------------------------------------------------------
void ComputeDamage::update(const int id_i, const int i) {
  const vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  const int jnum = PDconnections.size();
  const double maxConnections = jnum;
  if (maxConnections <= 0) {
    (*m_data)(i, m_indexDamage) = 1.;
    return;
  }

  double tot = 0;
  for (int jj = 0; jj < jnum; jj++) {
    const auto &con = PDconnections[jj];
    if (con.second[m_indexConnected] <= 0.5)
      continue;
    tot += 1;
  }

  (*m_data)(i, m_indexDamage) = 1. - tot / maxConnections;
  //    (*m_data)(i, m_indexDamage) = tot;
}
//------------------------------------------------------------------------------
void ComputeDamage::init(const int id_i, const int i) {
  (*m_data)(i, m_indexMaxPdConnections) =
      m_particles.pdConnections(id_i).size();
}
//------------------------------------------------------------------------------
}
