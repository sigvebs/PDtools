#include "force.h"

#include <stdio.h>

namespace PDtools
{
//------------------------------------------------------------------------------
Force::Force(PD_Particles &particles): m_particles(particles)
{

}
//------------------------------------------------------------------------------
Force::~Force()
{

}
//------------------------------------------------------------------------------
void Force::initialize(double E, double nu, double delta, int dim, double h)
{
    (void) E;
    (void) nu;
    (void) dim;
    (void) h;
    (void) delta;
}
//------------------------------------------------------------------------------
double Force::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    (void) idCol;
    return 0;
}
//------------------------------------------------------------------------------
void Force::calculatePotentialEnergy(const std::pair<int, int> &idCol,
                                     int indexPotential)
{
    (void) idCol;
    (void) indexPotential;
//    std::cerr << "ERROR: potential energy not implemented for this force" << std::endl;
}
//------------------------------------------------------------------------------
void Force::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    (void) idCol;
    (void) indexStress;
//    std::cerr << "ERROR: stress not implemented for this force" << std::endl;
}
//------------------------------------------------------------------------------
void Force::updateState()
{

}
//------------------------------------------------------------------------------
double Force::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    (void) idCol;
    (void) dt;

    return 0.;
}
//------------------------------------------------------------------------------
void Force::numericalInitialization(bool ni)
{
    m_numericalInitialization = ni;
}
//------------------------------------------------------------------------------
}
