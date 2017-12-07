#ifndef MESHTOPDPARTICES_H
#define MESHTOPDPARTICES_H

#include "PDtools/Mesh/pdmesh.h"

namespace PDtools {
class Grid;
class PD_Particles;
//------------------------------------------------------------------------------
PD_Particles convertMshToPdParticles(int dim, int interpolationDegree,
                                     const PdMesh &mesh, Grid &grid);
double areaTriangle(const mat &vertices);
}
#endif // MESHTOPDPARTICES_H
