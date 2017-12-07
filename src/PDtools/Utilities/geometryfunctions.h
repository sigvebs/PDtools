#ifndef GEOMETRYFUNCTIONS_H
#define GEOMETRYFUNCTIONS_H

#include "config.h"

class PD_Particles;
class PD_quadElement;

namespace PDtools {
//------------------------------------------------------------------------------
vec3 centroidOfQuad(PD_Particles &nodes, const PD_quadElement &quadElement);
//------------------------------------------------------------------------------
}
#endif // GEOMETRYFUNCTIONS_H
