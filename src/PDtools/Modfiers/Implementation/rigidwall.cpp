#include "rigidwall.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
RigidWall::RigidWall(Grid &grid, int orientationAxis, int topOrBottom,
                     double lc)
    : m_grid(grid), m_orientationAxis(orientationAxis),
      m_plane(topOrBottom), // 0 = bottom, 1 = top
      m_lc(lc) {
  switch (m_orientationAxis) {
  case 0:
    otherAxis = {1, 2};
    break;
  case 1:
    otherAxis = {0, 2};
    break;
  case 2:
    otherAxis = {0, 1};
    break;
  default:
    cerr << "Rigid Wall:: orientation:" << m_orientationAxis << " is not valid"
         << endl;
    exit(EXIT_FAILURE);
  }

  if (topOrBottom == 0)
    m_bottom = 1;
  if (topOrBottom == 1)
    m_top = 1;
}
//------------------------------------------------------------------------------
void RigidWall::initialize() {
  // Finding the gridpoints the beolngs to the wall
  // Assumes that the gridpoints are constant during the simulation
  const int me = m_grid.myRank();
  unordered_map<int, GridPoint> &gridpoints = m_grid.gridpoints();
  const arma::ivec3 nGrid = m_grid.nGrid();
  const int nx = nGrid(0);
  const int ny = nGrid(1);

  int topOrBottom = 0;
  if (m_plane == 1)
    topOrBottom = nGrid(m_orientationAxis) - 1;

  for (int i = 0; i < nGrid(otherAxis[0]); i++) {
    for (int j = 0; j < nGrid(otherAxis[1]); j++) {
      int x, y, z;
      switch (m_orientationAxis) {
      case 0:
        x = topOrBottom;
        y = i;
        z = j;
        break;
      case 1:
        x = i;
        y = topOrBottom;
        z = j;
        break;
      case 2:
        x = i;
        y = j;
        z = topOrBottom;
        break;
      default:
        exit(EXIT_FAILURE);
      }
      const int gridId = x + nx * y + nx * ny * z;
      const int belongsTo = m_grid.belongsTo(gridId);

      if (belongsTo == me) {
        m_gridpoints.push_back(&gridpoints.at(gridId));
        cout << "gridId:" << gridId << endl;
      }
    }
  }
}
//------------------------------------------------------------------------------
void RigidWall::evaluateStepOne() {
  m_boundary = m_grid.boundary();
  mat &R = m_particles->r();

  if (m_bottom) {
    const double boundary = m_boundary[m_orientationAxis].first - 0.48 * m_lc;
    for (GridPoint *gridPoint : m_gridpoints) {
      for (const pair<int, int> &idCol_i : gridPoint->particles()) {
        const int i = idCol_i.second;
        double x = R(i, m_orientationAxis);
        if (x < boundary) {
          R(i, m_orientationAxis) -= 2 * (x - boundary);
        }
      }
    }
  }
  if (m_top) {
    const double boundary = m_boundary[m_orientationAxis].second + 0.48 * m_lc;
    ;
    for (GridPoint *gridPoint : m_gridpoints) {
      for (const pair<int, int> &idCol_i : gridPoint->particles()) {
        const int i = idCol_i.second;
        double x = R(i, m_orientationAxis);
        if (x > boundary) {
          R(i, m_orientationAxis) -= 2 * (x - boundary);
        }
      }
    }
  }
}
//------------------------------------------------------------------------------
}
