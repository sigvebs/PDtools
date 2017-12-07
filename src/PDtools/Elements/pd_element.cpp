#include "pd_element.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Mesh/pdmesh.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Utilities/gaussianquadrature.h"
#include "PDtools/SavePdData/savepddata.h"

#define DEBUG_PRINT_ELEMENT_INIT 0
//------------------------------------------------------------------------------
namespace PDtools {
//------------------------------------------------------------------------------
PD_Particles initializeElementPd(const PdMesh &msh, const Grid &grid,
                                 const size_t quadratureDegree) {
  // For a review of quadrature rules for triangles and quads:
  // http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF

  const int dim = grid.dim();
  const int nTriangles = msh.nTriangles();
  const int nQuads = msh.nQuads();
  //    const int nElements = msh.nElements();
  //    const int nVertices = msh.nElements();

  const vector<array<double, 3>> &vertices = msh.vertices();
  const vector<vector<int>> &elements = msh.elements();
  const unordered_map<int, int> &IdToCol_vertices = msh.getIdToCol_vertices();
  const unordered_map<int, int> &IdToCol_elements = msh.getIdToCol_elements();

  //--------------------------------------------------------------------------
  // Initializing the Pdnodes
  //--------------------------------------------------------------------------
  const int nodesPerQuad = 4;
  const int nodesPerTri = 3;
  int nParticles = nodesPerQuad * nQuads + nodesPerTri * nTriangles;

  PD_Particles particles;
  particles.dim(dim);
  particles.maxParticles(nParticles);
  particles.nParticles(nParticles);
  particles.totParticles(nParticles);
  particles.initializeMatrices();
  particles.initializeElements(nTriangles, nQuads, quadratureDegree);

  unordered_map<int, int> &idToCol = particles.idToCol();
  arma::ivec &get_id = particles.colToId();
  arma::mat &r0 = particles.r0();
  arma::mat &r = particles.r();
  //    arma::mat & v = particles.v();
  //    arma::mat & data = particles.data();

  int elementId = 0;
  int pd_nodeId = 0;
  int pd_col = 0;
  //--------------------------------------------------------------------------
  // Setting the quadrature
  //--------------------------------------------------------------------------
  //    GaussLegendreQuad GaussLegendre_quadBasis(dim, quadratureDegree);
  orderedGrid GaussLegendre_quadBasis(dim, quadratureDegree);
  const int nQuadPoints = GaussLegendre_quadBasis.nPoints();

  //    mat gaussianPoints = GaussLegendre_quadBasis.gaussianPoints_2d();
  vec gaussianWeights = GaussLegendre_quadBasis.gaussianWeights_2d();
  mat shapeFunction = GaussLegendre_quadBasis.shapeFunction_2d();
  mat quadraturePoints = zeros(nQuadPoints, dim);
  vector<vector<double>> jakobi_1 = GaussLegendre_quadBasis.jakobi_1();
  vector<vector<double>> jakobi_2 = GaussLegendre_quadBasis.jakobi_2();

  particles.setShapeFunction(shapeFunction);
  //--------------------------------------------------------------------------
  double height = 1.;
  if (dim == 2) {
    const vector<pair<double, double>> &bondunary = grid.boundary();
    height = bondunary[2].second - bondunary[2].first;
  }

  switch (dim) {
  case 2:
    for (const vector<int> element : elements) {
      int counter = 0;

      switch (element.size()) {
      case 3: // TRIANGLE
      {
        array<size_t, 3> verticeIds;

        for (int vId : element) {
          int vCol = IdToCol_vertices.at(vId);
          const array<double, 3> &vert_i = vertices[vCol];

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
        double quad_x[4];
        double quad_y[4];
        int j = 0;
        for (int vId : element) {
          int vCol = IdToCol_vertices.at(vId);
          const array<double, 3> &vert_i = vertices[vCol];

          idToCol[pd_nodeId] = pd_col;
          get_id[pd_col] = pd_col;

          r(pd_col, 0) = vert_i[0];
          r(pd_col, 1) = vert_i[1];
          r(pd_col, 2) = vert_i[2];

          r0(pd_col, 0) = vert_i[0];
          r0(pd_col, 1) = vert_i[1];
          r0(pd_col, 2) = vert_i[2];

          quad_x[j] = vert_i[0];
          quad_y[j] = vert_i[1];

          verticeIds[counter] = pd_nodeId;
          pd_nodeId++;
          pd_col++;
          counter++;
          j++;
        }
//------------------------------------------------------------------
#if DEBUG_PRINT_ELEMENT_INIT
        double area = 0;
        double area_span = 0;
        for (int i = 0; i < 3; i++) {
          area_span += (quad_x[i] * quad_y[i + 1] - quad_x[i + 1] * quad_y[i]);
        }
        area_span += (quad_x[3] * quad_y[0] - quad_x[0] * quad_y[3]);
        area_span *= 0.5;
#endif
        // Setting the quadrature points
        vec weights(nQuadPoints);

        for (int i = 0; i < nQuadPoints; i++) {
          const vector<double> &dnda = jakobi_1[i];
          const vector<double> &dndb = jakobi_2[i];

          double dxda = 0;
          double dyda = 0;
          double dxdb = 0;
          double dydb = 0;

          double x = 0;
          double y = 0;
          for (int d = 0; d < 4; d++) {
            x += quad_x[d] * shapeFunction(i, d);
            y += quad_y[d] * shapeFunction(i, d);
            dxda += quad_x[d] * dnda[d];
            dxdb += quad_x[d] * dndb[d];
            dyda += quad_y[d] * dnda[d];
            dydb += quad_y[d] * dndb[d];
          }

          quadraturePoints(i, 0) = x;
          quadraturePoints(i, 1) = y;
          const double jakobiDeterminant = dxda * dydb - dxdb * dyda;
          weights(i) = height * gaussianWeights(i) * jakobiDeterminant;
#if DEBUG_PRINT_ELEMENT_INIT
          area += gaussianWeights(i) * jakobiDeterminant;
#endif
        }
#if DEBUG_PRINT_ELEMENT_INIT
        cout << shapeFunction << endl;
        cout.precision(8);
        quadraturePoints.col(0).raw_print(cout, "x = ");
        quadraturePoints.col(1).raw_print(cout, "y = ");

        cout << area << " " << area_span << endl;
#endif
        PD_quadElement qe(elementId, verticeIds);

        qe.setGuassianQuadraturePoints_initial(quadraturePoints);
        qe.setGuassianQuadraturePoints(quadraturePoints);
        qe.setGuassianQuadratureWeights(weights);
        particles.addQuad(qe);
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
  if (dim >= 1)
    saveparam_scale.push_back(pair<std::string, double>("x", 1.));
  if (dim >= 2)
    saveparam_scale.push_back(pair<std::string, double>("y", 1.));
  if (dim >= 3)
    saveparam_scale.push_back(pair<std::string, double>("z", 1.));

  saveparam_scale.push_back(pair<std::string, double>("volume", 1.));

  string saveParticlesPath = "geometry.xyz";
//  SaveParticles *saveParticles = new SaveParticles("xyz", saveparam_scale, false);
//  saveParticles->writeToFile(particles, saveParticlesPath);

  return particles;
}
//------------------------------------------------------------------------------
void updateElementQuadrature(PD_Particles &nodes) {
  // TODO: elements are not MPI ready - only the appropriate elements should be
  // updated
  const mat &R = nodes.r();
  const unordered_map<int, int> &idToCol = nodes.idToCol();
  vector<PD_quadElement> &quadElements = nodes.getQuadElements();
  const mat &shapeFunction = nodes.getShapeFunction();

  // Quad elements
  double quad_x[4];
  double quad_y[4];

  const size_t nQuadPoints = shapeFunction.n_rows;

  for (PD_quadElement &quadElement : quadElements) {
    const array<size_t, 4> &verticeIds = quadElement.verticeIds();
    int j = 0;

    mat &quadraturePoints = quadElement.guassianQuadraturePoints();

    for (size_t vId : verticeIds) {
      const int i = idToCol.at(vId);
      quad_x[j] = R(i, 0);
      quad_y[j] = R(i, 1);
      j++;
    }

    for (size_t i = 0; i < nQuadPoints; i++) {

      double x = 0;
      double y = 0;
      for (int d = 0; d < 4; d++) {
        x += quad_x[d] * shapeFunction(i, d);
        y += quad_y[d] * shapeFunction(i, d);
      }
      quadraturePoints(i, 0) = x;
      quadraturePoints(i, 1) = y;
    }
  }
}
//------------------------------------------------------------------------------
void printElement(const int element_id, PD_Particles &nodes) {
  (void)element_id;
  (void)nodes;
  /*
const mat &R = nodes.r();

const vector<PD_quadElement> &quadElements = nodes.getQuadElements();
const unordered_map<int, int> &idToElement = nodes.getIdToElement();

PD_quadElement quadElement = quadElements[idToElement.at(element_id)];
mat &quadraturePoints = quadElement.guassianQuadraturePoints();

const array<size_t, 4> &verticeIds = quadElement.verticeIds();

    int j = 0;
    for(size_t vId:verticeIds) {
        const int i = idToCol.at(vId);
        quad_x[j] = R(i,0);
        quad_y[j] = R(i,1);
        j++;
    }
    */
}
//------------------------------------------------------------------------------
}
