#ifndef LOADMESH_H
#define LOADMESH_H

#include "pdmesh.h"

//------------------------------------------------------------------------------
namespace PDtools {

enum MSH_ElementTypes {
  LINE = 1,
  TRIANGLE = 2,
  QUADRANGLE = 3,
  TETRAHEDRON = 4,
  HEXAHEDRON = 5,
  PRISM = 6,
  PYRAMID = 7
};

PdMesh loadMesh2d(string path);
PdMesh load_msh2d(string loadPath);
}
//------------------------------------------------------------------------------

#endif // LOADMESH_H
