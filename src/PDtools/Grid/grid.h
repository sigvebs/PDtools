#ifndef GRID_H
#define GRID_H

#include "config.h"

namespace PDtools {
class Particles;
class PD_Particles;

//------------------------------------------------------------------------------
class GridPoint {
public:
  GridPoint();
  GridPoint(int id, vec3 center, bool ghost);
  int id() const { return m_id; }
  const vec3 &center() const { return m_center; }
  bool isGhost() const { return m_ghost; }

  void setParticles(vector<pair<int, int>> particles) {
    m_particles = particles;
  }

  void addParticles(const vector<pair<int, int>> &particles) {
    m_particles.reserve(m_particles.size() + particles.size());
    m_particles.insert(m_particles.end(), particles.begin(), particles.end());
  }

  void addParticle(pair<int, int> p) { m_particles.push_back(p); }

  void addElement(const array<size_t, 2> &id_col);

  void removeParticle(const pair<int, int> &p) {
    m_particles.erase(std::remove(m_particles.begin(), m_particles.end(), p),
                      m_particles.end());
  }
  const vector<pair<int, int>> &particles() const { return m_particles; }
  void setNeighbours(vector<GridPoint *> neighbours) {
    m_neighbours = neighbours;
  }
  const vector<GridPoint *> &neighbours() const { return m_neighbours; }
  void clearParticles() { m_particles.clear(); }
  int ownedBy() const { return m_ownedBy; }
  void ownedBy(int ob) { m_ownedBy = ob; }

  void setNeighbourRanks(const vector<int> &neighbourRanks) {
    m_neighbourRanks = neighbourRanks;
  }

  const vector<int> &neighbourRanks() const { return m_neighbourRanks; }

  void setnGridId(const vector<int> &ids);

  const vector<int> &nGridId() const;

  vector<double> periodicShift() const;
  void setPeriodicShift(double shift, int d);

  int periodicNeighbourRank() const;

  void setPeriodicNeighbourRank(const int periodicNeighbourRank);

  vector<array<size_t, 2>> elements() const;

private:
  int m_id;
  vector<int> m_nGridId;
  int m_ownedBy = 0;
  vec3 m_center;
  bool m_ghost;
  vector<pair<int, int>> m_particles;
  vector<array<size_t, 2>> m_elements;
  vector<GridPoint *> m_neighbours;
  vector<int> m_neighbourRanks;
  int m_periodicNeighbourRank;
  vector<double> m_periodicShift = {0, 0, 0};
};
//------------------------------------------------------------------------------
class Grid {
private:
  int m_dim;
  double m_gridspacing;
  int m_myRank = 0;
  int m_nCores = 1;
  int m_nGrid[M_DIM];
  arma::ivec3 m_nGridArma;
  arma::ivec6 m_nGrid_with_boundary;
  arma::vec3 m_gridSpacing;
  std::unordered_map<int, GridPoint> m_gridpoints;
  std::vector<int> m_myGridPoints;
  std::vector<int> m_ghostGridIds;
  std::vector<int> m_periodicSendGridIds;
  std::vector<int> m_periodicReceiveGridIds;
  std::vector<int> m_neighbouringCores;
  std::vector<double> m_boundaryLength;
  std::vector<pair<double, double>> m_boundary;
  std::vector<std::array<double, 2>> m_boundary2;
  std::vector<pair<double, double>> m_originalBoundary;
  arma::ivec3 m_periodicBoundaries = {0, 0, 0};
  int m_counter = 0;
  vector<int> m_nCpuGrid;

  vector<int> m_boundaryGridPoints;
  enum m_enumCoordinates { X, Y, Z };

  double m_L0 = 1.0; // Scaling when loading particles

public:
  Grid();
  Grid(const vector<pair<double, double>> &boundary, double gridSpacing);
  Grid(const vector<pair<double, double>> &boundary, double gridSpacing,
       const ivec3 &periodicBoundaries);

  ~Grid();
  void initialize();
  void createGrid();
  void setNeighbours();
  int gridId(const vec3 &r) const;
  int gridId(const double (&r)[M_DIM]) const;
  int gridIdN(const vector<int> &nId) const;
  int particlesBelongsTo(const vec3 &r) const;
  void update();
  void placeParticlesInGrid(Particles &particles);
  void placeElementsInGrid(PD_Particles &nodes);
  void clearParticles();
  void clearElements();
  void clearAllParticles();
  void clearGhostParticles();
  void setIdAndCores(int myRank, int nCores);
  void setMyGridpoints();
  int belongsTo(const int gId) const;
  const vector<int> &myGridPoints() const { return m_myGridPoints; }

  unordered_map<int, GridPoint> &gridpoints() { return m_gridpoints; }

  void setOwnership();
  int myRank() { return m_myRank; }
  int nCores() { return m_nCores; }
  vector<int> &boundaryGridPoints();
  vector<int> &neighbouringCores();
  const std::vector<int> &ghostGrid() const;
  void setInitialPositionScaling(const double L0);
  double initialPositionScaling() const;
  const arma::ivec3 &nGrid() const;
  const vector<pair<double, double>> &boundary() const;
  vector<int> nCpuGrid() const;
  void setBoundaryGrid();
  std::vector<int> periodicSendGridIds() const;
  std::vector<int> periodicReceiveGridIds() const;
  std::vector<pair<double, double>> originalBoundary() const;

  int dim() const;
  void dim(int dim);
};
//------------------------------------------------------------------------------
// Inline functions
//------------------------------------------------------------------------------
inline vector<int> &Grid::boundaryGridPoints() { return m_boundaryGridPoints; }
inline vector<int> &Grid::neighbouringCores() { return m_neighbouringCores; }
//------------------------------------------------------------------------------
// Other grid dependent functions
//------------------------------------------------------------------------------
void updateVerletList(const std::string &verletId, Particles &particles,
                      Grid &grid, double radius);

//------------------------------------------------------------------------------
}
#endif // GRID_H
