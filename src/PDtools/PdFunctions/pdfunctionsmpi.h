#ifndef PDFUNCTIONSMPI_H
#define PDFUNCTIONSMPI_H

//------------------------------------------------------------------------------
namespace PDtools {

class PD_Particles;
class Grid;
class Force;
class Modifier;

void exchangeGhostParticles(Grid &grid, PD_Particles &particles);
void exchangeGhostParticles_boundary(Grid &grid, PD_Particles &particles);
void exchangeInitialGhostParticles(Grid &grid, PD_Particles &particles);
void exchangeInitialGhostParticles_boundary(Grid &grid,
                                            PD_Particles &particles);
void updateGrid(Grid &grid, PD_Particles &particles, const bool ADR = false);
void updateModifierLists(Modifier &modifier, PD_Particles &particles,
                         int counter);
void exchangeInitialPeriodicBoundaryParticlesInitial(Grid &grid,
                                                     PD_Particles &particles,
                                                     bool sendPDdata);
void exchangePeriodicBoundaryParticlesInitial(Grid &grid,
                                              PD_Particles &particles,
                                              bool sendPDdata);
}
//------------------------------------------------------------------------------
#endif // PDFUNCTIONSMPI_H
