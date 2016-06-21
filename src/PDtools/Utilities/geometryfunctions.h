#ifndef GEOMETRYFUNCTIONS_H
#define GEOMETRYFUNCTIONS_H

#include "Particles/particles.h"
#include "Particles/pd_particles.h"
#include <armadillo>

using namespace std;
using namespace arma;

namespace PDtools
{
//------------------------------------------------------------------------------


vec3 centroidOfQuad(PD_Particles &nodes, const PD_quadElement &quadElement);

//------------------------------------------------------------------------------
}
#endif // GEOMETRYFUNCTIONS_H
