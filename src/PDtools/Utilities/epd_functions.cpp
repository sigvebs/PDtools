#include "epd_functions.h"

#include "PDtools/Elements/pd_element.h"
#include "PDtools/Grid/grid.h"
#include "Particles/pd_particles.h"
#include <fstream>

namespace PDtools {
//------------------------------------------------------------------------------
int saveElementQuadrature(PD_Particles &nodes, Grid &grid,
                          const std::string &savePath) {
  std::ofstream outStream;
  outStream.open(savePath.c_str());

  const vector<pair<double, double>> &boundary = grid.originalBoundary();
  const double scaling = 1.0;
  const double xMin = scaling * boundary[0].first;
  const double xMax = scaling * boundary[0].second;
  const double yMin = scaling * boundary[1].first;
  const double yMax = scaling * boundary[1].second;
  const double zMin = scaling * boundary[2].first;
  const double zMax = scaling * boundary[2].second;

  const mat &shapeFunction = nodes.getShapeFunction();
  const size_t nQuadPoints = shapeFunction.n_rows;
  vector<PD_quadElement> &quadElements = nodes.getQuadElements();
  const int nQuadElements = nQuadPoints * quadElements.size();

  outStream << "ITEM: TIMESTEP" << endl;
  outStream << 0 << endl;
  outStream << "ITEM: NUMBER OF ATOMS" << endl;
  outStream << nQuadElements << endl;
  outStream << "ITEM: BOX BOUNDS pp pp pp" << endl;
  outStream << xMin << ' ' << xMax << "\n";
  outStream << yMin << ' ' << yMax << "\n";
  outStream << zMin << ' ' << zMax << "\n";

  outStream << "ITEM: ATOMS";
  outStream << " id x y qId" << endl;

  outStream.setf(std::ios::scientific);
  outStream.precision(14);

  int id_i = 0;
  for (PD_quadElement &quadElement : quadElements) {
    const int pId = quadElement.id();
    mat &quadraturePoints = quadElement.guassianQuadraturePoints();

    for (size_t i = 0; i < nQuadPoints; i++) {
      outStream << id_i << " " << quadraturePoints(i, 0) << " "
                << quadraturePoints(i, 1) << " " << pId << endl;
      id_i++;
    }
  }
  outStream.close();

  return 1;
}
//------------------------------------------------------------------------------
}
