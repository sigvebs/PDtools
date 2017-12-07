#include "geometryfunctions.h"

#include "PDtools/Elements/pd_element.h"
#include "PDtools/Particles/particles.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
vec3 centroidOfQuad(PD_Particles &nodes, const PD_quadElement &quadElement) {
  const mat &R = nodes.r();
  const ivec &idToCol = nodes.getIdToCol_v();

  const array<size_t, 4> vertexIds = quadElement.verticeIds();
  vec3 r = {0, 0, 0};
  for (size_t vertex_id : vertexIds) {
    const int i = idToCol[vertex_id];
    r += R.row(i).t();
  }

  return r / vertexIds.size();
}
//------------------------------------------------------------------------------
}
