#ifndef FORCE_H
#define FORCE_H

#include <unordered_map>
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
    int m_dim = 3;
    double m_E;
    double m_nu;
    double m_h;
    double m_delta;
    double m_lc;

    arma::mat & m_r;
    arma::mat & m_r0;
    arma::mat & m_F;
    arma::mat & m_data;
    std::unordered_map<int, int> & m_pIds;

    enum enum_coordinates{X, Y, Z};
    static constexpr double THRESHOLD = 2.2204e-016;
public:
    Force(PD_Particles &particles);

    virtual ~Force();

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    virtual void
    calculateForces(const std::pair<int, int> & idCol ) = 0;

    virtual double
    calculatePotentialEnergyDensity(const std::pair<int, int> & idCol);

    virtual void
    calculatePotentialEnergy(const std::pair<int, int> & idCol,
                                          int indexPotential);
    virtual double
    calculateBondEnergy(const std::pair<int, int> & idCol,
                                          std::pair<int, std::vector<double>> & con);

    virtual void
    calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);
    virtual void
    updateState();

    virtual void
    updateState(const std::pair<int, int> & idCol);

    virtual double
    calculateStableMass(const std::pair<int, int> & idCol,
                                     double dt);
    void
    numericalInitialization(bool ni);

    void
    setDim(int dim);

    virtual void
    applySurfaceCorrection(double strain = 0.001);

    virtual void
    applyStrainCorrection(double strain);

    virtual void
    applyShearCorrection(double shear=0.001);
};
//------------------------------------------------------------------------------
}
#endif // FORCE_H
