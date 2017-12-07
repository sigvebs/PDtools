#include "computemaxstretch.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
ComputeMaxStretch::ComputeMaxStretch(PD_Particles &particles)
    : ComputeProperty(particles) {
  m_indexStretch = m_particles.getPdParamId("stretch");
  m_indexMaxStretch = m_particles.registerParameter("max_stretch");
  m_data = &m_particles.data();
}
//------------------------------------------------------------------------------
ComputeMaxStretch::~ComputeMaxStretch() {}
//------------------------------------------------------------------------------
void ComputeMaxStretch::update(const int id_i, const int i) {
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id_i);
  double sMax = 0;
  for (auto &con : PDconnections) {
    double s = con.second[m_indexStretch];
    sMax = std::max(s, sMax);
    //        sMax = max(fabs(s), sMax);
  }
  (*m_data)(i, m_indexMaxStretch) = sMax;
}
//------------------------------------------------------------------------------
void ComputeMaxStretch::init(const int id_i, const int i) {
  (void)id_i;
  (void)i;
}
//------------------------------------------------------------------------------
}
