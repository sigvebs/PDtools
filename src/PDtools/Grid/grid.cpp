#include "grid.h"

#include "Domain/domain.h"
#include "Particles/particles.h"
#include "Particles/pd_particles.h"
#include <map>
namespace PDtools
{
//------------------------------------------------------------------------------
Grid::Grid()
{

}
//------------------------------------------------------------------------------
Grid::Grid(const Domain &domain, double gridspacing):
    m_dim(domain.dim), m_gridspacing(gridspacing),
    m_periodicBoundaries(domain.periodicBoundaries())
{
    m_boundaryLength = domain.boundaryLength();
    m_boundary = domain.boundaries();
}
//------------------------------------------------------------------------------
Grid::Grid(const vector<pair<double, double>> &boundary, double gridspacing):
    m_dim(boundary.size()),
    m_gridspacing(gridspacing),
    m_boundary(boundary)
{
    for(int d=0; d<m_dim; d++)
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
    if(m_dim <= 1)
    {
        ny = 1;
        y_len = m_gridspacing;
        m_gridSpacing(1) = 1.;
    }
    if(m_dim <= 2)
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
            }
        }
    }

    // Boundary conditions
    for(int d=0; d<3; d++)
    {
        ivec3 xyz = {0, 0, 0};

        for(int point_0:boundary_points[d])
        {
            vector<int> inn_p = {0, 1, 2};
            inn_p.erase(inn_p.begin() + d);

            const ivec p1 = join_cols(inner_points[inn_p[0]], boundary_points[inn_p[0]]);
            for(int point_1:p1)
            {
                const ivec p2 = join_cols(inner_points[inn_p[1]], boundary_points[inn_p[1]]);

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

    setOwnership();
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
        for(int d=0; d<m_dim; d++)
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

    for(int d=0; d<DIM; d++)
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

#if DIM == 2
    return i(X) + m_nGrid(X)*i(Y);
#else
    return i(X) + m_nGrid(X)*i(Y) + m_nGrid(X)*m_nGrid(Y)*i(Z);
#endif
}
//------------------------------------------------------------------------------
int Grid::particlesBelongsTo(const vec3 &r) const
{
    const int gId = gridId(r);
//    return m_gridpoints.at(gId)->ownedBy();
    return belongsTo(gId);
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
    const ivec & colToId = particles.colToId();

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        const int id = colToId.at(i);
        const vec3 &r = R.row(i).t();
        const int gId = gridId(r);
        const pair<int, int> id_pos(id, i);
#ifdef USE_OPENMP
#pragma omp critical
#endif
        m_gridpoints[gId]->addParticle(id_pos);
    }
}
//------------------------------------------------------------------------------
void Grid::clearParticles()
{
    for(int id:m_myGridPoints)
    {
        GridPoint* gp = m_gridpoints[id];
        gp->clearParticles();
    }
    clearGhostParticles();
}
//------------------------------------------------------------------------------
void Grid::clearAllParticles()
{
    for(pair<int, GridPoint*> id_gridpoint:m_gridpoints)
    {
        id_gridpoint.second->clearParticles();
    }
}
//------------------------------------------------------------------------------
void Grid::clearGhostParticles()
{
    for(int id:m_ghostGridIds)
    {
        GridPoint* gp = m_gridpoints[id];
        gp->clearParticles();
    }
}
//------------------------------------------------------------------------------
void Grid::setIdAndCores(int myRank, int nCores)
{
    m_myRank = myRank;
    m_nCores = nCores;
}
//------------------------------------------------------------------------------
void Grid::setMyGridpoints()
{
    vector<int> boundaryGridPoints;
    vector<int> ghostGridIds;
    vector<int> neighbouringCores;

    for(auto &gridPoint:m_gridpoints)
    {
        if(gridPoint.second->ownedBy() == m_myRank)
        {
            const int id = gridPoint.first;
            m_myGridPoints.push_back(id);

            // Setting boundary cells and cpu-ids
            const vector<GridPoint*> & neighbours = gridPoint.second->neighbours();
            vector<int> neighbourRanks;

            for(const GridPoint *neighbour:neighbours)
            {
                const int neighbourRank = neighbour->ownedBy();
                if(neighbourRank != m_myRank)
                {
                    neighbourRanks.push_back(neighbourRank);
                    ghostGridIds.push_back(neighbour->id());

                    bool found = false;
                    for(int i:neighbouringCores)
                    {
                        if(neighbourRank == i)
                            found = true;
                    }
                    if(!found)
                        neighbouringCores.push_back(neighbourRank);
                }
            }

            if(neighbourRanks.size() > 0)
            {
                // Only picking the uniqe ids
                sort( neighbourRanks.begin(), neighbourRanks.end() );
                neighbourRanks.erase( unique( neighbourRanks.begin(),
                                              neighbourRanks.end() ),
                                      neighbourRanks.end() );

                boundaryGridPoints.push_back(id);
                gridPoint.second->setNeighbourRanks(neighbourRanks);
            }
        }
    }
    sort( neighbouringCores.begin(), neighbouringCores.end() );
    m_boundaryGridPoints = boundaryGridPoints;
    m_ghostGridIds = ghostGridIds;
    m_neighbouringCores = neighbouringCores;
}
//------------------------------------------------------------------------------
int Grid::belongsTo(const int gId) const
{
    const double nPoints = m_gridpoints.size();
    return (m_nCores*(gId + 1.) - 1.)/nPoints;
}
//------------------------------------------------------------------------------
void Grid::setOwnership()
{
    for(const auto &gridPoint:m_gridpoints)
    {
        const int id = gridPoint.first;
        gridPoint.second->ownedBy(belongsTo(id));
    }
}
//------------------------------------------------------------------------------
const std::vector<int> &Grid::ghostGrid()
{
    return m_ghostGridIds;
}
//------------------------------------------------------------------------------
void Grid::setInitialPositionScaling(const double L0)
{
    m_L0 = L0;
}
//------------------------------------------------------------------------------
double Grid::initialPositionScaling()
{
    return m_L0;
}
//------------------------------------------------------------------------------
const arma::ivec3 &Grid::nGrid()
{
    return m_nGrid;
}
//------------------------------------------------------------------------------
const vector<pair<double, double> > &Grid::boundary()
{
    return m_boundary;
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

    const int verletId = particles.getVerletId(verletStringId);
    const unordered_map<int, GridPoint*> &gridpoints = grid.gridpoints();
    const mat & R = particles.r();
    const vector<int> &mygridPoints = grid.myGridPoints();
    const int dim = grid.dim();
    particles.clearVerletList(verletId);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0; i<mygridPoints.size(); i++)
    {
        double dr;
        const int gridId = mygridPoints.at(i);
        const GridPoint & gridPoint = *gridpoints.at(gridId);

        for(const pair<int, int> & idCol_i:gridPoint.particles())
        {
            const int id_i = idCol_i.first;
            const int i = idCol_i.second;
            vector<int> verletList;
            for(const pair<int, int> & idCol_j:gridPoint.particles())
            {
                const int id_j = idCol_j.first;
                const int j = idCol_j.second;

                if(id_i == id_j)
                    continue;

                double drSquared = 0;
                for(int d=0;d<dim; d++)
                {
                    dr = R(i, d) - R(j, d);
                    drSquared += dr*dr;
                }

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
                    const int id_j = idCol_j.first;
                    const int j = idCol_j.second;
                    double drSquared = 0;

                    for(int d=0;d<dim; d++)
                    {
                        dr = R(i, d) - R(j, d);
                        drSquared += dr*dr;
                    }

                    if(drSquared < radiusSquared)
                    {
                        verletList.push_back(id_j);
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
