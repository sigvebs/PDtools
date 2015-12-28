#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <math.h>
#include <unordered_map>
#include <armadillo>

namespace PDtools
{
using namespace std;
using namespace arma;

// Forward declerations
class Domain;
class Particles;

//------------------------------------------------------------------------------
class GridPoint
{
public:
    GridPoint();
    GridPoint(int id, vec3 center, bool ghost);

    int id()
    {
        return m_id;
    }

    const vec3 & center() const
    {
        return m_center;
    }

    bool isGhost() const
    {
        return m_ghost;
    }

    void setParticles(vector<pair<int, int>> particles)
    {
        m_particles = particles;
    }

    void addParticles(const vector<pair<int, int>> &particles)
    {
        m_particles.reserve(m_particles.size() + particles.size());
        m_particles.insert(m_particles.end(), particles.begin(), particles.end());
    }

    void addParticle(pair<int, int> p)
    {
        m_particles.push_back(p);
    }

    const vector<pair<int, int>> & particles() const
    {
        return m_particles;
    }

    void setNeighbours(vector<GridPoint*> neighbours)
    {
        m_neighbours = neighbours;
    }

    const vector<GridPoint*> & neighbours() const
    {
        return m_neighbours;
    }

    void clearParticles()
    {
        m_particles.clear();
    }

protected:
    int m_id;
    vec3 m_center;
    bool m_ghost;
    vector<pair<int,int>> m_particles;
    vector<GridPoint*> m_neighbours;
    vector<int> m_n;
};
//------------------------------------------------------------------------------
class Grid
{
public:
    int dim;

protected:
    double m_gridspacing;
    arma::ivec3 m_nGrid;
    arma::ivec6 m_nGrid_with_boundary;
    arma::vec3 m_gridSpacing;
    std::unordered_map<int, GridPoint*> m_gridpoints;
    std::vector<int> m_myGridPoints;
    std::vector<double> m_boundaryLength;
    std::vector<pair<double, double>> m_boundary;
    arma::ivec3 m_periodicBoundaries = {0, 0, 0};
    int m_counter = 0;

    enum m_enumCoordinates{X, Y, Z};

public:
    Grid();
    Grid(const Domain & domain, double gridspacing);
    Grid(const vector<pair<double, double> > &boundary, double gridSpacing);
    Grid(const vector<pair<double, double> > &boundary, double gridSpacing,
                  const ivec3 &periodicBoundaries);

    ~Grid();
    void initialize();
    void createGrid();
    void setNeighbours();
    int gridId(const vec3 & r) const;
    void update();
    void placeParticlesInGrid(Particles &particles);
    void clearParticles();

    const vector<int> &myGridPoints() const
    {
        return m_myGridPoints;
    }

    unordered_map<int, GridPoint*> & gridpoints()
    {
        return m_gridpoints;
    }
};
//------------------------------------------------------------------------------
// Other grid dependent functions
//------------------------------------------------------------------------------
void updateVerletList(const std::string &verletId,
                                  Particles & particles,
                                  Grid & grid,
                                  double radius);

//------------------------------------------------------------------------------
}
#endif // GRID_H
