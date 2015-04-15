#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <math.h>
#include <armadillo>
#include <unordered_map>
#include "Domain/domain.h"
#include "Particles/particles.h"

using namespace std;
using namespace arma;
namespace PDtools
{
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

    void addParticle(pair<int, int> p)
    {
        m_particles.push_back(p);
    }

    const vector<pair<int, int>> & particles() const
    {
        return m_particles;
    }

    void setNeighbours(vector<pair<int, GridPoint&>> neighbours)
    {
        m_neighbours = neighbours;
    }

    const vector<pair<int, GridPoint&>> & neighbours() const
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
    vector<pair<int, GridPoint&>> m_neighbours;
    vector<int> m_n;
};
//------------------------------------------------------------------------------
class Grid
{
public:
    const int dim;

protected:
    double m_gridspacing;
    ivec3 m_nGrid;
    ivec6 m_nGrid_with_boundary;
    vec3 m_gridSpacing;
    unordered_map<int, GridPoint> m_gridpoints;
    vector<double> m_boundaryLength;
    vector<pair<double, double>> m_boundary;
    ivec3 m_periodicBoundaries = {0, 0, 0};

    enum m_enumCoordinates{X, Y, Z};

public:
    Grid(const Domain & domain, double gridspacing);
    Grid(const vector<pair<double, double> > &boundary, double gridSpacing);
    Grid(const vector<pair<double, double> > &boundary, double gridSpacing,
                  const ivec3 &periodicBoundaries);

//    ~Grid();
    void initialize();
    void createGrid();
    void setNeighbours();
    int gridId(const vec3 & r) const;
    void update();
    void placeParticlesInGrid(Particles &particles);
    void clearParticles();

    unordered_map<int, GridPoint> & gridpoints()
    {
        return m_gridpoints;
    }
};
//------------------------------------------------------------------------------
}
#endif // GRID_H
