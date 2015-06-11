#ifndef FORCE_H
#define FORCE_H

#include <armadillo>

namespace PDtools
{
// Forward declerations
class PD_Particles;

//------------------------------------------------------------------------------
class Force
{
protected:
    PD_Particles &m_particles;
    bool m_numericalInitialization = false;

public:
    Force(PD_Particles &particles);

    virtual ~Force();

    virtual void initialize(double E, double nu, double delta, int dim, double h);

    virtual void calculateForces(const std::pair<int, int> & idCol ) = 0;


    virtual double calculatePotentialEnergyDensity(const std::pair<int, int> & idCol);
    virtual void calculatePotentialEnergy(const std::pair<int, int> & idCol,
                                          int indexPotential);

    virtual void calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);
    virtual void updateState();

    virtual double calculateStableMass(const std::pair<int, int> & idCol,
                                     double dt);

    static constexpr double THRESHOLD = 2.2204e-016;

    enum enum_coordinates{X, Y, Z};

    void numericalInitialization(bool ni);

};
//------------------------------------------------------------------------------
}
#endif // FORCE_H
