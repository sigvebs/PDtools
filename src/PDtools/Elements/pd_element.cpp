#include "pd_element.h"

#include "Particles/pd_particles.h"
#include "SavePdData/savepddata.h"
#include "Grid/grid.h"
#include "Mesh/pdmesh.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
//PD_elementsAndNodes::PD_elementsAndNodes(Grid &grid):
//    m_dim(grid.dim()), m_grid(grid)
//{

//}
////------------------------------------------------------------------------------
//void PD_elementsAndNodes::initialize(const PdMesh &msh)
//{

//}

//------------------------------------------------------------------------------
PD_Particles initializeElementPd(const PdMesh &msh, const Grid &grid, const size_t quadratureDegree)
{
    const int dim = grid.dim();
    const int nTriangles = msh.nTriangles();
    const int nQuads = msh.nQuads();
    const int nElements = msh.nElements();
    const int nVertices = msh.nElements();

    const vector<array<double, 3> > &vertices = msh.vertices();
    const vector<vector<int> > &elements = msh.elements();
    const unordered_map<int, int> & IdToCol_vertices = msh.getIdToCol_vertices();
    const unordered_map<int, int> & IdToCol_elements = msh.getIdToCol_elements();

    //--------------------------------------------------------------------------
    // Initializing the Pdnodes
    //--------------------------------------------------------------------------
    const int nodesPerQuad = 4;
    const int nodesPerTri = 3;
    int nParticles = nodesPerQuad*nQuads + nodesPerTri*nTriangles;

    PD_Particles particles;
    particles.dim(dim);
    particles.maxParticles(nParticles);
    particles.nParticles(nParticles);
    particles.totParticles(nParticles);
    particles.initializeMatrices();
    particles.initializeElements(nTriangles, nQuads, quadratureDegree);

    unordered_map<int, int> & idToCol = particles.idToCol();
    arma::ivec & get_id = particles.colToId();
    arma::mat & r0 = particles.r0();
    arma::mat & r = particles.r();
    arma::mat & v = particles.v();
    arma::mat & data = particles.data();

    int elementId = 0;
    int pd_nodeId = 0;
    int pd_col = 0;

    //--------------------------------------------------------------------------
    // Setting the quadrature
    //--------------------------------------------------------------------------
    GaussLegendreQuad GaussLegendre_quadBasis(dim, quadratureDegree);
    const int nQuadPoints = GaussLegendre_quadBasis.nPoints();

    mat gaussianPoints = GaussLegendre_quadBasis.gaussianPoints_2d();
    vec gaussianWeights = GaussLegendre_quadBasis.gaussianWeights_2d();
    mat shapeFunction = GaussLegendre_quadBasis.shapeFunction_2d();
    mat quadraturePoints = zeros(nQuadPoints, dim);
    //--------------------------------------------------------------------------

    switch (dim) {
    case 2:
    for(const vector<int> element:elements) {
        int counter = 0;

        switch (element.size()) {
        case 3: // TRIANGLE
        {
            array<size_t, 3> verticeIds;

            for(int vId:element) {
                int vCol = IdToCol_vertices.at(vId);
                const array<double, 3> & vert_i = vertices[vCol];

                idToCol[pd_nodeId] = pd_col;
                get_id[pd_col] = pd_col;

                r(pd_col, 0) = vert_i[0];
                r(pd_col, 1) = vert_i[1];
                r(pd_col, 2) = vert_i[2];

                r0(pd_col, 0) = vert_i[0];
                r0(pd_col, 1) = vert_i[1];
                r0(pd_col, 2) = vert_i[2];

                verticeIds[counter] = pd_nodeId;
                pd_nodeId++;
                pd_col++;
                counter++;
            }
            particles.addTri(PD_triElement(elementId, verticeIds));
            elementId++;
            break;
        }
        case 4: // QUAD
        {
            array<size_t, 4> verticeIds;

            for(int vId:element) {
                int vCol = IdToCol_vertices.at(vId);
                const array<double, 3> & vert_i = vertices[vCol];

                idToCol[pd_nodeId] = pd_col;
                get_id[pd_col] = pd_col;

                r(pd_col, 0) = vert_i[0];
                r(pd_col, 1) = vert_i[1];
                r(pd_col, 2) = vert_i[2];

                r0(pd_col, 0) = vert_i[0];
                r0(pd_col, 1) = vert_i[1];
                r0(pd_col, 2) = vert_i[2];

                verticeIds[counter] = pd_nodeId;
                pd_nodeId++;
                pd_col++;
                counter++;
            }


//            // Setting the quadrature points
//            for(int i=0; i<nQuadPoints; i++) {
//                const double x =
//                quadraturePoints(i, 0) = ;
//            }

            particles.addQuad(PD_quadElement(elementId, verticeIds));
            elementId++;
            break;
        }
        default:
            break;
        }
    }
        break;
    default:
        break;
    }

    particles.nParticles(pd_col);
    cout << pd_col << " " << nParticles << endl;

    //--------------------------------------------------------------------------
    // Saving the data
    //--------------------------------------------------------------------------
    vector<pair<string, double>> saveparam_scale;
    saveparam_scale.push_back(pair<std::string, double>("id", 1.));
    if(dim >= 1)
        saveparam_scale.push_back(pair<std::string, double>("x", 1.));
    if(dim >= 2)
        saveparam_scale.push_back(pair<std::string, double>("y", 1.));
    if(dim >= 3)
        saveparam_scale.push_back(pair<std::string, double>("z", 1.));

    saveparam_scale.push_back(pair<std::string, double>("volume", 1.));

    string saveParticlesPath = "geometry.xyz";
    SaveParticles *saveParticles = new SaveParticles("xyz", saveparam_scale, false);
    saveParticles->writeToFile(particles, saveParticlesPath);

    return particles;

}

//------------------------------------------------------------------------------
//template<size_t T>
//const array<size_t, T> &PD_element::verticeIds() const
//{
//    return m_verticeIds;
//}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
}
