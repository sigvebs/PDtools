#include "computegridid.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Grid/grid.h"

namespace PDtools {
//------------------------------------------------------------------------------
ComputeGridId::ComputeGridId(PD_Particles &particles, Grid &grid)
    : ComputeProperty(particles), m_grid(grid), m_r(particles.r()),
      m_data(particles.data()) {
  m_indexgrid = m_particles.registerParameter("gridId");
}
//------------------------------------------------------------------------------
ComputeGridId::~ComputeGridId() {}
//------------------------------------------------------------------------------
void ComputeGridId::update(const int id_i, const int i) {
  (void)id_i;
  const vec3 &l_r = m_r.row(i).t();
  const int gridId = m_grid.gridId(l_r);
  m_data(i, m_indexgrid) = gridId;
}
//------------------------------------------------------------------------------
}
