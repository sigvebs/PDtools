#include "grid.h"

#include <array>
#include <map>
#include "PDtools/PdFunctions/pdfunctions.h"
#include "Particles/particles.h"
#include "Particles/pd_particles.h"
#include "Utilities/geometryfunctions.h"


namespace PDtools
{
//------------------------------------------------------------------------------
vector<int> Grid::nCpuGrid() const
{
    return m_nCpuGrid;
}
//------------------------------------------------------------------------------
std::vector<int> Grid::periodicSendGridIds() const
{
    return m_periodicSendGridIds;
}
//------------------------------------------------------------------------------
std::vector<int> Grid::periodicReceiveGridIds() const
{
    return m_periodicReceiveGridIds;
}
//------------------------------------------------------------------------------
std::vector<pair<double, double> > Grid::originalBoundary() const
{
    return m_originalBoundary;
}
//------------------------------------------------------------------------------
int Grid::dim() const
{
    return m_dim;
}
//------------------------------------------------------------------------------
void Grid::dim(int dim)
{
    m_dim = dim;
}
//------------------------------------------------------------------------------
Grid::Grid()
{

}
//------------------------------------------------------------------------------
Grid::Grid(const vector<pair<double, double>> &boundary, double gridspacing):
    m_dim(boundary.size()),
    m_gridspacing(gridspacing),
    m_boundary(boundary),
    m_originalBoundary(boundary)
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
    setOwnership();
    setNeighbours();
    setBoundaryGrid();
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
    if(m_dim <= 1) {
        ny = 1;
        y_len = m_gridspacing;
        m_gridSpacing(1) = 1.;
    }
    if(m_dim <= 2) {
        nz = 1;
        z_len = m_gridspacing;
        m_gridSpacing(2) = 1.;
    }

    m_nGrid[0] = nx;
    m_nGrid[1] = ny;
    m_nGrid[2] = nz;

    vector<ivec> inner_points;
    vector<ivec> boundary_points;

    for(int d=0;d<M_DIM; d++) {
        if(m_periodicBoundaries[d]) {
            inner_points.push_back(linspace<ivec>(1, m_nGrid[d], m_nGrid[d]));
            boundary_points.push_back({0, m_nGrid[d] + 1});
            m_nGrid[d] += 2;
            m_boundary[d].first -= m_gridSpacing(d);
            m_boundary[d].second += m_gridSpacing(d);
            m_boundaryLength[d] +=  2*m_gridSpacing(d);
        } else {
            inner_points.push_back(linspace<ivec>(0, m_nGrid[d]-1, m_nGrid[d]));
            boundary_points.push_back({});
        }
    }

    // Setting the grid
    for(int x:inner_points[X]) {
        for(int y:inner_points[Y]) {
            for(int z:inner_points[Z]) {
                // id = x + nx*y + nx*ny*z;
                // Center of gridpoint
                const vec3 center = {
                    m_boundary[X].first + m_gridSpacing(X)*(x + 0.5)
                    ,m_boundary[Y].first + m_gridSpacing(Y)*(y + 0.5)
                    ,m_boundary[Z].first + m_gridSpacing(Z)*(z + 0.5)};

                const int id = gridId(center);
                const vector<int> nId = {x, y, z};
                m_gridpoints[id] = GridPoint(id, center, false);
                m_gridpoints.at(id).setnGridId(nId);
            }
        }
    }

    // Boundary conditions
    for(int d=0; d<M_DIM; d++) {
        Col<int> xyz = {0, 0, 0};

        for(int point_0:boundary_points[d]) {
            vector<int> inn_p = {0, 1, 2};
            inn_p.erase(inn_p.begin() + d);

            const ivec p1 = join_cols(inner_points[inn_p[0]], boundary_points[inn_p[0]]);
            for(int point_1:p1) {
                const ivec p2 = join_cols(inner_points[inn_p[1]], boundary_points[inn_p[1]]);

                for(int point_2:p2) {
                    xyz(d) = point_0;
                    xyz(inn_p[0]) = point_1;
                    xyz(inn_p[1]) = point_2;

                    // Center of gridpoint
                    const vec3 center = {
                        m_boundary[X].first + m_gridSpacing(X)*(xyz(X) + 0.5)
                        ,m_boundary[Y].first + m_gridSpacing(Y)*(xyz(Y) + 0.5)
                        ,m_boundary[Z].first + m_gridSpacing(Z)*(xyz(Z) + 0.5)};

                    const int id = gridId(center);
                    const vector<int> nId = {xyz(X), xyz(Y), xyz(Z)};
                    m_gridpoints[id] = GridPoint(id, center, true);
                    m_gridpoints.at(id).setnGridId(nId);
                }
            }
        }
    }

    m_nGridArma = {0, 0, 0};
    for(int d=0;d<M_DIM; d++) {
        m_nGridArma(d) = m_nGrid[d];
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
    for(int a:shift) {
        for(int b:shift) {
            for(int c:shift) {
                permuations.push_back(ivec3({a, b, c}));
            }
        }
    }

    // Setting the neighbours --------------------------------------------------
    for(const pair<const int, GridPoint> &id_gridpoint:m_gridpoints)
    {
        const int id = id_gridpoint.first;
        GridPoint & gridpoint = m_gridpoints[id];

        const vec3 center = gridpoint.center();
        ivec3 l_gridId = {0, 0, 0};
        for(int d=0; d<m_dim; d++)
            l_gridId(d) = int((center(d) - m_boundary[d].first)/m_gridSpacing(d));

        vector<GridPoint*> neighbours;
        ivec3 gridId_neigh = {0, 0, 0};

        for(ivec3 & shift_xyz:permuations) {
            gridId_neigh = l_gridId + shift_xyz;
            vector<bool> outOfBounds = {false, false, false};

            for(int d=0; d<3; d++) {
                if(m_periodicBoundaries[d]) {
                    if(gridId_neigh(d) == m_nGrid[d]-1) {
                        gridId_neigh[d] = 0;
                        outOfBounds[d] = false;
                    }
                    else if(gridId_neigh[d] == 0) {
                        gridId_neigh[d] = m_nGrid[d]-1;
                        outOfBounds[d] = false;
                    }
                }
                else {
                    if(gridId_neigh(d) >= m_nGrid[d] || gridId_neigh[d] < 0)
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

            neighbours.push_back(&m_gridpoints[id_neighbour]);
        }
        gridpoint.setNeighbours(neighbours);
    }
}
//------------------------------------------------------------------------------
int Grid::gridId(const vec3 &r) const
{
    int i[M_DIM];

    for(int d=0; d<M_DIM; d++) {
        i[d] = int((r(d) - m_boundary[d].first)/m_gridSpacing(d));

        // Handling the boundary extremals
        if(i[d] >= m_nGrid[d]) {
            i[d]  = m_nGrid[d] - 1;
        }
        else if(i[d]  < 0) {
            i[d]  = 0;
        }
    }

#if M_DIM == 2
    return i[X] + m_nGrid[X]*i[Y];
#else
    return i[X] + m_nGrid[X]*i[Y] + m_nGrid[X]*m_nGrid[Y]*i[Z];
#endif
}
//------------------------------------------------------------------------------
int Grid::gridIdN(const vector<int> &xyz) const
{
    const vec3 center = {
        m_boundary[X].first + m_gridSpacing(X)*(xyz[X] + 0.5)
        ,m_boundary[Y].first + m_gridSpacing(Y)*(xyz[Y] + 0.5)
        ,m_boundary[Z].first + m_gridSpacing(Z)*(xyz[Z] + 0.5)};

    return gridId(center);
}
//------------------------------------------------------------------------------
int Grid::particlesBelongsTo(const vec3 &r) const
{
    const int gId = gridId(r);
    return m_gridpoints.at(gId).ownedBy();
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
        const size_t id = colToId.at(i);
        const vec3 &r = R.row(i).t();
        const size_t gId = gridId(r);
        const pair<int, int> id_pos(id, i);
#ifdef USE_OPENMP
#pragma omp critical
#endif
        m_gridpoints[gId].addParticle(id_pos);
    }
}
//------------------------------------------------------------------------------
void Grid::placeElementsInGrid(PD_Particles &nodes)
{
    clearElements();
    const vector<PD_quadElement> & quadElements = nodes.getQuadElements();

    // Quad elements
    for(size_t i=0; i< quadElements.size(); i++) {
        const size_t e_id = quadElements[i].id();
        const PD_quadElement & quadElement = quadElements[i];
        const vec3 &r = centroidOfQuad(nodes, quadElement);
        const size_t gId = gridId(r);
        m_gridpoints[gId].addElement({e_id, i});
    }
}
//------------------------------------------------------------------------------
void Grid::clearParticles()
{
    for(int id:m_myGridPoints)
    {
        GridPoint& gp = m_gridpoints[id];
        gp.clearParticles();
    }
    clearGhostParticles();
}
//------------------------------------------------------------------------------
void Grid::clearElements()
{
//    for(int id:m_myGridPoints)
//    {
//        GridPoint* gp = m_gridpoints[id];
//        gp->clearParticles();
//    }
//    clearGhostParticles();
}
//------------------------------------------------------------------------------
void Grid::clearAllParticles()
{
    for(pair<const int, GridPoint> &id_gridpoint:m_gridpoints) {
        id_gridpoint.second.clearParticles();
    }
}
//------------------------------------------------------------------------------
void Grid::clearGhostParticles()
{
    for(int id:m_ghostGridIds)
    {
        GridPoint& gp = m_gridpoints[id];
        gp.clearParticles();
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

    for(auto &gridPoint:m_gridpoints) {
        GridPoint & gp = gridPoint.second;

        if(gp.ownedBy() == m_myRank) {
            const int id = gridPoint.first;
            m_myGridPoints.push_back(id);

            // Setting boundary cells and cpu-ids
            const vector<GridPoint*> & neighbours = gp.neighbours();
            vector<int> neighbourRanks;

            for(const GridPoint *neighbour:neighbours) {
                const int neighbourRank = neighbour->ownedBy();
                if(neighbourRank != m_myRank) {
                    neighbourRanks.push_back(neighbourRank);
                    ghostGridIds.push_back(neighbour->id());

                    bool found = false;
                    for(int i:neighbouringCores) {
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
                gp.setNeighbourRanks(neighbourRanks);
            }
        }
    }
    sort( neighbouringCores.begin(), neighbouringCores.end() );
    sort(m_myGridPoints.begin(), m_myGridPoints.end());
    m_boundaryGridPoints = boundaryGridPoints;
    m_ghostGridIds = ghostGridIds;
    m_neighbouringCores = neighbouringCores;
}
//------------------------------------------------------------------------------
int Grid::belongsTo(const int gId) const
{
    return m_gridpoints.at(gId).ownedBy();
//    const double nPoints = m_gridpoints.size();
//    return (m_nCores*(gId + 1.) - 1.)/nPoints;
}
//------------------------------------------------------------------------------
void Grid::setOwnership()
{
    m_nCpuGrid = optimalConfigurationCores(m_nCores, m_boundaryLength, m_dim);

    int rank;
    for(auto &gridPoint:m_gridpoints) {
        GridPoint & gp = gridPoint.second;
        vector<int> n = gridPoint.second.nGridId();
#if USE_MPI
        int i[M_DIM];

        for(int d=0; d<M_DIM; d++) {
            i[d] = int(((n[d])/float(m_nGrid[d]))*(m_nCpuGrid[d]));
        }

#if M_DIM == 2
        rank = i[X] + m_nCpuGrid[X]*i[Y];
#else
        rank = i[X] + m_nCpuGrid[X]*i[Y] + m_nCpuGrid[X]*m_nCpuGrid[Y]*i[Z];
#endif
//        const double nPoints = m_gridpoints.size();
//        rank = (m_nCores*(id + 1.) - 1.)/nPoints;
        gp.ownedBy(rank);
#else
        gp.ownedBy(0);
#endif
    }
}
//------------------------------------------------------------------------------
void Grid::setBoundaryGrid()
{
    for(int d=0; d<m_dim; d++) {
        if(m_periodicBoundaries[d]) {
            vector<int> yz;
            switch (d) {
            case 0:
                yz = {1, 2};
                break;
            case 1:
                yz = {0, 2};
                break;
            case 2:
                yz = {0, 1};
                break;
            }

            double shift = m_boundaryLength[d] - 2*m_gridSpacing(d);

            int x0 = 1;
            int xn = m_nGrid[d] - 2;
            for(int y=0; y<m_nGrid[yz[0]]; y++) {
                for(int z=0; z<m_nGrid[yz[1]]; z++) {
                    // Left boundary
                    vector<int> leftIdN = {0, 0, 0};
                    leftIdN[d] = x0;
                    leftIdN[yz[0]] = y;
                    leftIdN[yz[1]] = z;
                    const int leftId = gridIdN(leftIdN);
                    GridPoint& leftGp = m_gridpoints[leftId];
                    const int leftOwnedBy = leftGp.ownedBy();
                    vector<int> leftIdN_0 = leftIdN;
                    leftIdN_0[d] = x0-1;
                    const int leftId_0 = gridIdN(leftIdN_0);
                    GridPoint& leftGp_0 = m_gridpoints[leftId_0];
                    const int leftOwnedBy_0 = leftGp_0.ownedBy();

                    // Right boundary
                    vector<int> rightIdN = {0, 0, 0};
                    rightIdN[d] = xn;
                    rightIdN[yz[0]] = y;
                    rightIdN[yz[1]] = z;
                    const int rightId = gridIdN(rightIdN);
                    GridPoint& rightGp = m_gridpoints[rightId];
                    const int rightOwnedBy = rightGp.ownedBy();
                    vector<int> rightIdN_0 = rightIdN;
                    rightIdN_0[d] = xn + 1;
                    const int rightId_0 = gridIdN(rightIdN_0);
                    GridPoint& rightGp_0 = m_gridpoints[rightId_0];
                    const int rightOwnedBy_0 = rightGp_0.ownedBy();

                    // Setting the shifts
                    leftGp.setPeriodicShift(shift, d);
                    rightGp.setPeriodicShift(-shift, d);
                    leftGp_0.setPeriodicShift(shift, d);
                    rightGp_0.setPeriodicShift(-shift, d);
                    leftGp.setPeriodicNeighbourRank(rightOwnedBy);
                    rightGp.setPeriodicNeighbourRank(leftOwnedBy);
                    leftGp_0.setPeriodicNeighbourRank(rightOwnedBy_0);
                    rightGp_0.setPeriodicNeighbourRank(leftOwnedBy_0);

                    if(leftOwnedBy == m_myRank) {
                        m_periodicSendGridIds.push_back(leftId);
                        m_periodicReceiveGridIds.push_back(leftId_0);
//                        cout << m_myRank << "->" << rightOwnedBy << " L left:" << leftId << " right:" << rightId << " "
//                             << x0 << " "<< y << " " << z << endl;
                    }
                    if(rightOwnedBy == m_myRank) {
                        m_periodicSendGridIds.push_back(rightId);
                        m_periodicReceiveGridIds.push_back(rightId_0);
//                        cout << m_myRank  << "->" << leftOwnedBy << " R left:" << leftId << " right:" << rightId << " "
//                             << xn << " "<< y << " " << z << endl;
                    }
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
const std::vector<int> &Grid::ghostGrid() const
{
    return m_ghostGridIds;
}
//------------------------------------------------------------------------------
void Grid::setInitialPositionScaling(const double L0)
{
    m_L0 = L0;
}
//------------------------------------------------------------------------------
double Grid::initialPositionScaling() const
{
    return m_L0;
}
//------------------------------------------------------------------------------
const arma::ivec3 &Grid::nGrid() const
{
    return m_nGridArma;
}
//------------------------------------------------------------------------------
const vector<pair<double, double> > &Grid::boundary() const
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
void GridPoint::addElement(const array<size_t, 2> &id_col)
{
    m_elements.push_back(id_col);
}
//------------------------------------------------------------------------------
void GridPoint::setnGridId(const vector<int> &ids)
{
    m_nGridId = ids;
}
//------------------------------------------------------------------------------
const vector<int> &GridPoint::nGridId() const
{
    return m_nGridId;
}
//------------------------------------------------------------------------------
vector<double> GridPoint::periodicShift() const
{
    return m_periodicShift;
}
//------------------------------------------------------------------------------
void GridPoint::setPeriodicShift(double shift, int d)
{
    m_periodicShift[d] = shift;
}
//------------------------------------------------------------------------------
int GridPoint::periodicNeighbourRank() const
{
    return m_periodicNeighbourRank;
}
//------------------------------------------------------------------------------
void GridPoint::setPeriodicNeighbourRank(const int periodicNeighbourRank)
{
    m_periodicNeighbourRank = periodicNeighbourRank;
}
//------------------------------------------------------------------------------
vector<array<size_t, 2> > GridPoint::elements() const
{
    return m_elements;
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
    const unordered_map<int, GridPoint> &gridpoints = grid.gridpoints();
    const mat & R = particles.r();
    const vector<int> &mygridPoints = grid.myGridPoints();
    const int dim = grid.dim();
    particles.clearVerletList(verletId);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0; i<mygridPoints.size(); i++) {
        double dr;
        const int gridId = mygridPoints.at(i);
        const GridPoint & gridPoint = gridpoints.at(gridId);
        if(gridPoint.isGhost())
            continue;

        for(const pair<int, int> & idCol_i:gridPoint.particles()) {
            const int id_i = idCol_i.first;
            const int i = idCol_i.second;
            vector<int> verletList;
            for(const pair<int, int> & idCol_j:gridPoint.particles()) {
                const int id_j = idCol_j.first;
                const int j = idCol_j.second;

                if(id_i == id_j)
                    continue;

                double drSquared = 0;
                for(int d=0;d<dim; d++) {
                    dr = R(i, d) - R(j, d);
                    drSquared += dr*dr;
                }

                if(drSquared < radiusSquared) {
                    verletList.push_back(id_j);
                }
            }
            // Neighbouring cells

            const vector<GridPoint*> & neighbours = gridPoint.neighbours();

            for(const GridPoint *neighbour:neighbours) {
                for(const pair<int, int> & idCol_j:neighbour->particles()) {
                    const int id_j = idCol_j.first;
                    const int j = idCol_j.second;
                    double drSquared = 0;

                    for(int d=0;d<dim; d++) {
                        dr = R(i, d) - R(j, d);
                        drSquared += dr*dr;
                    }

                    if(drSquared < radiusSquared) {
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
