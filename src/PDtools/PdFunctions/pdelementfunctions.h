#ifndef PDELEMENTFUNCTIONS_H
#define PDELEMENTFUNCTIONS_H

namespace PDtools {
class PD_Particles;
class Grid;
//------------------------------------------------------------------------------

void setPdElementConnections(PD_Particles &discretization, Grid &grid,
                             const double delta);
//------------------------------------------------------------------------------
}

#endif // PDELEMENTFUNCTIONS_H
