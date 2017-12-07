#include "pdelementfunctions.h"

#include "Grid/grid.h"
#include "PDtools/Elements/pd_element.h"
#include "Particles/pd_particles.h"
#include "Utilities/boostgeometry_settings.h"
#include "config.h"

#define DEBUG_PRINT_ELEMENTS 0

namespace PDtools {
//------------------------------------------------------------------------------
void setPdElementConnections(PD_Particles &discretization, Grid &grid,
                             const double delta) {
  const unordered_map<int, GridPoint> &gridpoints = grid.gridpoints();
  const vector<int> &mygridPoints = grid.myGridPoints();
  const mat &R = discretization.r();
  const ivec &idToCol = discretization.getIdToCol_v();
  const vector<PD_quadElement> &quadElements = discretization.getQuadElements();

  // The order is important!
  discretization.registerPdParameter("connected");
  discretization.registerPdParameter("overlap");
  vector<size_t> pCols;
  pCols.reserve(4);

  const polygon_2d baseCicle = createCircle(point_2d(0, 0), delta);
//  double A0 = M_PI * delta * delta;

  const double delta2 = delta * delta;

#ifdef USE_OPENMP
#pragma omp parallel for private(pCols)
#endif
  for (size_t i = 0; i < mygridPoints.size(); i++) {
    double dx, dy;
    const size_t gridId = mygridPoints.at(i);
    const GridPoint &gridPoint = gridpoints.at(gridId);

    for (const pair<int, int> &idCol_i : gridPoint.particles()) {
      int id_i = idCol_i.first;
      int col_i = idCol_i.second;
      unordered_map<int, vector<double>> connections;
      vector<pair<int, vector<double>>> connectionsVector;

      int totoalNumberPolygons = 0;
      int intersectingPolygons = 0;

      for (const array<size_t, 2> &idCol_j : gridPoint.elements()) {
        const int j = idCol_j[1];
        const PD_quadElement &quadElement = quadElements[j];
        const array<size_t, 4> &pIds = quadElement.verticeIds();

        bool isInside = false;

        for (int k = 0; k < 4; k++) {
          const int col_j = idToCol[pIds[k]];
          pCols[k] = col_j;
          if (isInside)
            continue;
          dx = R(col_i, 0) - R(col_j, 0);
          dx *= dx;
          dy = R(col_i, 1) - R(col_j, 1);
          dy *= dy;
          const double dr2 = dx + dy;
          if (dr2 < delta2)
            isInside = true;
        }
        if (isInside == false)
          continue;

        const int element_id = idCol_j[0];
        const double overlap = 1;

        vector<double> connectionData;
        connectionData.push_back(1.0);     // Connected
        connectionData.push_back(overlap); // Overlap factor
        connections[element_id] = connectionData;
        connectionsVector.push_back(
            pair<int, vector<double>>(element_id, connectionData));
        intersectingPolygons++;
        totoalNumberPolygons++;
      }

      // Neighbouring cells
      const vector<GridPoint *> &neighbours = gridPoint.neighbours();
      for (const GridPoint *neighbour : neighbours) {
        for (const array<size_t, 2> &idCol_j : neighbour->elements()) {
          const int j = idCol_j[1];
          const PD_quadElement &quadElement = quadElements[j];
          const array<size_t, 4> &pIds = quadElement.verticeIds();

          bool isInside = false;

          for (int k = 0; k < 4; k++) {
            const int col_j = idToCol[pIds[k]];
            pCols[k] = col_j;
            if (isInside)
              continue;
            dx = R(col_i, 0) - R(col_j, 0);
            dx *= dx;
            dy = R(col_i, 1) - R(col_j, 1);
            dy *= dy;
            const double dr2 = dx + dy;
            if (dr2 < delta2)
              isInside = true;
          }
          if (isInside == false)
            continue;

          const int element_id = idCol_j[0];
          const double overlap = 1;

          vector<double> connectionData;
          connectionData.push_back(1.0);     // Connected
          connectionData.push_back(overlap); // Overlap factor
          connections[element_id] = connectionData;
          connectionsVector.push_back(
              pair<int, vector<double>>(element_id, connectionData));
          intersectingPolygons++;

          /*
          polygon_2d polygon_j;
          polygon_j.outer().push_back(point_2d(R(pCols[0], 0), R(pCols[0], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[3], 0), R(pCols[3], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[2], 0), R(pCols[2], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[1], 0), R(pCols[1], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[0], 0), R(pCols[0], 1)));

          std::deque<polygon_2d> output;
          boost::geometry::intersection(interactionPolygon_i, polygon_j,
          output);
          if(output.size() > 0) {
              const int element_id = idCol_j[0];
              const polygon_2d & intersection = output[0];
              const double area_intersection = bg::area(intersection);
              const double area_p = bg::area(polygon_j);
              const double overlap = area_intersection/area_p;

              vector<double> connectionData;
              connectionData.push_back(1.0); // Connected
              connectionData.push_back(overlap); // Overlap factor
              connections[element_id] = connectionData;
              connectionsVector.push_back(pair<int, vector<double>>(element_id,
          connectionData) );
              intersectingPolygons++;
              tot_area += area_intersection;
          }
          */
          totoalNumberPolygons++;
        }
      }

#ifdef USE_OPENMP
#pragma omp critical
      { particles.setPdConnections(id_i, connectionsVector); }
#else
      discretization.setPdConnections(id_i, connectionsVector);
#endif
#if DEBUG_PRINT_ELEMENTS
      if (tot_area / A0 < 0.98) {

        totoalNumberPolygons = 0;
        intersectingPolygons = 0;
        double tot_area = 0;
        for (const array<size_t, 2> &idCol_j : gridPoint.elements()) {
          const int j = idCol_j[1];
          const PD_quadElement &quadElement = quadElements[j];
          const array<size_t, 4> &pIds = quadElement.verticeIds();

          bool isInside = false;

          for (int k = 0; k < 4; k++) {
            const int col_j = idToCol[pIds[k]];
            pCols[k] = col_j;
            if (isInside)
              continue;
            dx = R(col_i, 0) - R(col_j, 0);
            dx *= dx;
            dy = R(col_i, 1) - R(col_j, 1);
            dy *= dy;
            const double dr2 = dx + dy;
            if (dr2 < delta2)
              isInside = true;
          }
          //                    if(isInside == false)
          //                        continue;

          polygon_2d polygon_j;
          polygon_j.outer().push_back(point_2d(R(pCols[0], 0), R(pCols[0], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[3], 0), R(pCols[3], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[2], 0), R(pCols[2], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[1], 0), R(pCols[1], 1)));
          polygon_j.outer().push_back(point_2d(R(pCols[0], 0), R(pCols[0], 1)));

          //                    bool intersects_ij =
          //                    bg::intersects(interactionPolygon_i, polygon_j);

          //                    if (!intersects_ij)
          //                        continue;

          std::deque<polygon_2d> output;
          boost::geometry::intersection(interactionPolygon_i, polygon_j,
                                        output);
          if (output.size() > 0) {
            const polygon_2d &intersection = output[0];
            const double area_intersection = bg::area(intersection);

            intersectingPolygons++;
            tot_area += area_intersection;
            std::cout << bg::wkt<polygon_2d>(intersection) << std::endl;
          }
          std::cout << bg::wkt<polygon_2d>(polygon_j) << std::endl;

          totoalNumberPolygons++;
        }

        // Neighbouring cells
        const vector<GridPoint *> &neighbours = gridPoint.neighbours();
        for (const GridPoint *neighbour : neighbours) {
          for (const array<size_t, 2> &idCol_j : neighbour->elements()) {
            const int j = idCol_j[1];
            const PD_quadElement &quadElement = quadElements[j];
            const array<size_t, 4> &pIds = quadElement.verticeIds();

            bool isInside = false;

            for (int k = 0; k < 4; k++) {
              const int col_j = idToCol[pIds[k]];
              pCols[k] = col_j;
              if (isInside)
                continue;
              dx = R(col_i, 0) - R(col_j, 0);
              dx *= dx;
              dy = R(col_i, 1) - R(col_j, 1);
              dy *= dy;
              const double dr2 = dx + dy;
              if (dr2 < delta2)
                isInside = true;
            }
            //                        if(isInside == false)
            //                            continue;

            polygon_2d polygon_j;
            polygon_j.outer().push_back(
                point_2d(R(pCols[0], 0), R(pCols[0], 1)));
            polygon_j.outer().push_back(
                point_2d(R(pCols[3], 0), R(pCols[3], 1)));
            polygon_j.outer().push_back(
                point_2d(R(pCols[2], 0), R(pCols[2], 1)));
            polygon_j.outer().push_back(
                point_2d(R(pCols[1], 0), R(pCols[1], 1)));
            polygon_j.outer().push_back(
                point_2d(R(pCols[0], 0), R(pCols[0], 1)));

            //                        bool intersects_ij =
            //                        bg::intersects(interactionPolygon_i,
            //                        polygon_j);
            //                        if (!intersects_ij)
            //                            continue;

            std::deque<polygon_2d> output;
            boost::geometry::intersection(interactionPolygon_i, polygon_j,
                                          output);
            if (output.size() > 0) {
              const polygon_2d &intersection = output[0];
              const double area_intersection = bg::area(intersection);
              intersectingPolygons++;
              tot_area += area_intersection;

              std::cout << bg::wkt<polygon_2d>(intersection) << std::endl;
            }
            std::cout << bg::wkt<polygon_2d>(polygon_j) << std::endl;
            totoalNumberPolygons++;
          }
        }
        //                bg::correct(p);
        //                std::cout << std::boolalpha << bg::is_valid(p) <<
        //                "\n";
        cout << double(intersectingPolygons) / double(totoalNumberPolygons)
             << "\t " << tot_area / A0 << endl;
        cout << "----------------" << endl;
        //                    bool intersects_ij =
        //                    bg::intersects(interactionPolygon_i, polygon_j);
        //                    if (!intersects_ij)
        //                        continue;
      }
#endif
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
}
