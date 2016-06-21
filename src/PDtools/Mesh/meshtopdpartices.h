#ifndef MESHTOPDPARTICES_H
#define MESHTOPDPARTICES_H

#include "Mesh/pdmesh.h"
#include "Particles/pd_particles.h"
#include "Utilities/gaussianquadrature.h"

namespace PDtools
{
class Grid;

//------------------------------------------------------------------------------
PD_Particles convertMshToPdParticles(int dim, int interpolationDegree, const PdMesh &mesh, Grid &grid);
double areaTriangle(const mat &vertices);

}
#endif // MESHTOPDPARTICES_H
