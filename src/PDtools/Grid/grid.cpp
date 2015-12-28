#include "grid.h"

#include "Domain/domain.h"
#include "Particles/particles.h"
#include <map>

namespace PDtools
{
//------------------------------------------------------------------------------
Grid::Grid()
{

}
//------------------------------------------------------------------------------
Grid::Grid(const Domain &domain, double gridspacing):
    dim(domain.dim), m_gridspacing(gridspacing),
    m_periodicBoundaries(domain.periodicBoundaries())
{
    m_boundaryLength = domain.boundaryLength();
    m_boundary = domain.boundaries();
}
//------------------------------------------------------------------------------
Grid::Grid(const vector<pair<double, double>> &boundary, double gridspacing):
    dim(boundary.size()),
    m_gridspacing(gridspacing),
    m_boundary(boundary)
{
    for(int d=0; d<dim; d++)
    {
        m_boundaryLength.push_back(m_boundary[d].second - m_boundary[d].first);
    }
}
//------------------------------------------------------------------------------
Grid::Grid(const vector<pair<double, double> > &boundary, double gridspacing,
           const ivec3 &periodicBoundaries):
    Grid(boundary, gridspacing)
{
    m_periodicBoundaries = periodicBoundaries;
}
//------------------------------------------------------------------------------
void Grid::initialize()
{
    update();
    createGrid();
    setNeighbours();
}
//------------------------------------------------------------------------------
void Grid::createGrid()
{
    double x_len = m_boundaryLength[0];
    double y_len = m_boundaryLength[1];
    double z_len = m_boundaryLength[2];

    // Finding the optimal length
    int nx = floor(x_len/m_gridspacing) > 0 ? floor(x_len/m_gridspacing) : 1;
    int ny = floor(y_len/m_gridspacing) > 0 ? floor(y_len/m_gridspacing) : 1;
    int nz = floor(z_len/m_gridspacing) > 0 ? floor(z_len/m_gridspacing) : 1;


    m_gridSpacing = {x_len/nx, y_len/ny, z_len/nz };
    if(dim <= 1)
    {
        ny = 1;
        y_len = m_gridspacing;
        m_gridSpacing(1) = 1.;
    }
    if(dim <= 2)
    {
        nz = 1;
        z_len = m_gridspacing;
        m_gridSpacing(2) = 1.;
    }

    m_nGrid = {nx, ny, nz};
    vector<ivec> inner_points;
    vector<ivec> boundary_points;

    for(int d=0;d<3; d++)
    {
        if(m_periodicBoundaries[d])
        {
            inner_points.push_back(linspace<ivec>(1, m_nGrid(d), m_nGrid(d)));
            boundary_points.push_back({0, m_nGrid(d) + 1});
            m_nGrid(d) += 2;
            m_boundary[d].first -= m_gridSpacing(d);
            m_boundary[d].second += m_gridSpacing(d);
        }
        else
        {
            inner_points.push_back(linspace<ivec>(0, m_nGrid(d)-1, m_nGrid(d)));
            boundary_points.push_back({});
        }
    }

    // Setting the grid
    for(int x:inner_points[X])
    {
        for(int y:inner_points[Y])
        {
            for(int z:inner_points[Z])
            {
                // id = x + nx*y + nx*ny*z;
                // Center of gridpoint
                const vec3 center = {
                    m_boundary[X].first + m_gridSpacing(X)*(x + 0.5)
                    ,m_boundary[Y].first + m_gridSpacing(Y)*(y + 0.5)
                    ,m_boundary[Z].first + m_gridSpacing(Z)*(z + 0.5)};

                const int id = gridId(center);
                m_gridpoints[id] = new GridPoint(id, center, false);
                m_myGridPoints.push_back(id);
            }
        }
    }
    cout << m_myGridPoints.size() << endl;

    // Boundary conditions
    for(int d=0; d<3; d++)
    {
        ivec3 xyz = {0, 0, 0};

        for(int point_0:boundary_points[d])
        {
            vector<int> inn_p = {0, 1, 2};
            inn_p.erase(inn_p.begin() + d);

            ivec p1 = join_cols(inner_points[inn_p[0]], boundary_points[inn_p[0]]);
            for(int point_1:p1)
            {
                ivec p2 = join_cols(inner_points[inn_p[1]], boundary_points[inn_p[1]]);

                for(int point_2:p2)
                {
                    xyz(d) = point_0;
                    xyz(inn_p[0]) = point_1;
                    xyz(inn_p[1]) = point_2;

                    // Center of gridpoint
                    const vec3 center = {
                        m_boundary[X].first + m_gridSpacing(X)*(xyz(X) + 0.5)
                        ,m_boundary[Y].first + m_gridSpacing(Y)*(xyz(Y) + 0.5)
                        ,m_boundary[Z].first + m_gridSpacing(Z)*(xyz(Z) + 0.5)};

                    const int id = gridId(center);
                    m_gridpoints[id] = new GridPoint(id, center, true);
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
void Grid::setNeighbours()
{
    vector<int> shift;
    shift.push_back(-1);
    shift.push_back(0);
    shift.push_back(1);

    // Creating all possible permuations of (-1,0,1)^3
    vector<ivec3> permuations;
    for(int a:shift)
    {
        for(int b:shift)
        {
            for(int c:shift)
            {
                permuations.push_back(ivec3({a, b, c}));
            }
        }
    }

    // Setting the neighbours --------------------------------------------------
    for(pair<int, GridPoint*> id_gridpoint:m_gridpoints)
    {
        const int id = id_gridpoint.first;
        GridPoint & gridpoint = *m_gridpoints[id];

        if(gridpoint.isGhost())
            continue;

        const vec3 center = gridpoint.center();
        ivec3 l_gridId = {0, 0, 0};
        for(int d=0; d<dim; d++)
            l_gridId(d) = int((center(d) - m_boundary[d].first)/m_gridSpacing(d));

        vector<GridPoint*> neighbours;
        ivec3 gridId_neigh = {0, 0, 0};


        for(ivec3 & shift_xyz:permuations)
        {
            gridId_neigh = l_gridId + shift_xyz;
            vector<bool> outOfBounds = {false, false, false};

            for(int d=0; d<3; d++)
            {
                if(m_periodicBoundaries[d])
                {
                    if(gridId_neigh(d) == m_nGrid(d)-1)
                    {
                        gridId_neigh[d] = 0;
                        outOfBounds[d] = false;
                    }
                    else if(gridId_neigh[d] == 0)
                    {
                        gridId_neigh[d] = m_nGrid(d)-1;
                        outOfBounds[d] = false;
                    }
                }
                else
                {
                    if(gridId_neigh(d) >= m_nGrid(d) || gridId_neigh[d] < 0)
                        outOfBounds[d] = true;
                }
            }

            if(outOfBounds[X] || outOfBounds[Y] || outOfBounds[Z])
                continue;

            vec3 coord_neigh;
            coord_neigh = center + shift_xyz % m_gridSpacing;
            int id_neighbour = gridId(coord_neigh);

            if(id_neighbour == id)
                continue;

            neighbours.push_back(m_gridpoints[id_neighbour]);
        }
        gridpoint.setNeighbours(neighbours);
    }
}
//------------------------------------------------------------------------------
int Grid::gridId(const vec3 &r) const
{
    ivec3 i;

    for(int d=0; d<3; d++)
//        for(int d=0; d<dim; d++)
    {
        i(d) = int((r(d) - m_boundary[d].first)/m_gridSpacing(d));

        // Handling the boundary extremals
        if(i(d) >= m_nGrid(d))
        {
            i(d) = m_nGrid(d) - 1;
        }
        else if(i(d) < 0)
        {
            i(d) = 0;
        }
    }
    return i(X) + m_nGrid(X)*i(Y) + m_nGrid(X)*m_nGrid(Y)*i(Z);
}
//------------------------------------------------------------------------------
void Grid::update()
{
    //    m_boundaryLength = m_domain.boundaryLength();
    //    m_boundary = m_domain.boundaries();
}

//------------------------------------------------------------------------------
void Grid::placeParticlesInGrid(Particles &particles)
{
    clearParticles();
    const mat & R = particles.r();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> id_pos(i, i);
        const vec3 &r = R.row(id_pos.second).t();
        int gId = gridId(r);
#ifdef USE_OPENMP
#pragma omp critical
#endif
        m_gridpoints[gId]->addParticle(id_pos);
    }
}
//------------------------------------------------------------------------------
void Grid::clearParticles()
{
    for(pair<int, GridPoint*> id_gridpoint:m_gridpoints)
    {
        id_gridpoint.second->clearParticles();
    }
}
//------------------------------------------------------------------------------
Grid::~Grid()
{
    //    m_gridpoints.clear();
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Gridpoint functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
GridPoint::GridPoint()
{
}
//------------------------------------------------------------------------------
GridPoint::GridPoint(int id, vec3 center, bool ghost):
    m_id(id), m_center(center), m_ghost(ghost)
{
}
//------------------------------------------------------------------------------
// Other grid dependent functions
//------------------------------------------------------------------------------

void updateVerletList(const string &verletStringId,
                      Particles & particles,
                      Grid & grid,
                      double radius)
{
    double radiusSquared = radius*radius;

    int verletId = particles.getVerletId(verletStringId);
    const unordered_map<int, GridPoint*> &gridpoints = grid.gridpoints();
    const mat & R = particles.r();
    const vector<int> &mygridPoints = grid.myGridPoints();
    particles.clearVerletList(verletId);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0; i<mygridPoints.size(); i++)
    {
        double dx, dy, dz;
        int gridId = mygridPoints.at(i);
        const GridPoint & gridPoint = *gridpoints.at(gridId);

        for(const pair<int, int> & idCol_i:gridPoint.particles())
        {
            int id_i = idCol_i.first;
            vector<int> verletList;
            const vec & r_i = R.row(idCol_i.second);

            for(const pair<int, int> & idCol_j:gridPoint.particles())
            {
                int id_j = idCol_j.first;
                if(id_i == id_j)
                    continue;

                const vec & r_j = R.row(idCol_j.second);
                dx = r_i(0) - r_j(0);
                dy = r_i(1) - r_j(1);
                dz = r_i(2) - r_j(2);

                double drSquared = dx*dx + dy*dy + dz*dz;

                if(drSquared < radiusSquared)
                {
                    verletList.push_back(id_j);
                }
            }

            // Neighbouring cells
            const vector<GridPoint*> & neighbours = gridPoint.neighbours();

            for(const GridPoint *neighbour:neighbours)
            {
                for(const pair<int, int> & idCol_j:neighbour->particles())
                {
                    const vec & r_j = R.row(idCol_j.second);
                    dx = r_i(0) - r_j(0);
                    dy = r_i(1) - r_j(1);
                    dz = r_i(2) - r_j(2);

                    double drSquared = dx*dx + dy*dy + dz*dz;

                    if(drSquared < radiusSquared)
                    {
                        verletList.push_back(idCol_j.first);
                    }
                }
            }

#ifdef USE_OPENMP
#pragma omp critical
            {
                particles.setVerletList(id_i, verletList, verletId);
            }
#else
            particles.setVerletList(id_i, verletList, verletId);
#endif
        }
    }
}

//------------------------------------------------------------------------------
}
