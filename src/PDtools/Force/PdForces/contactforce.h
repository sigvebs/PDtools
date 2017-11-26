#ifndef CONTACTFORCE_H
#define CONTACTFORCE_H

#include <unordered_map>
#include "PDtools/Force/force.h"

namespace PDtools
{
class Grid;
//------------------------------------------------------------------------------
class ContactForce : public Force
{
private:
    Grid & m_grid;
    int m_steps;
    double m_spacing;
    double m_scaling;
    double m_verletRadius;
    double m_forceScaling;
    double m_scaling_dr0;
    double m_velocityScaling;
    int m_indexVolume;
    int m_verletUpdateFrq = 10;
    int m_verletListId;
    int m_indexMicromodulus;
    int m_indexRadius;
    int m_iDr0;
    int m_indexConnected;
    int m_dim;
    std::unordered_map<int, int> *m_idToCol ;
public:
    ContactForce(PD_Particles &particles, Grid &grid, double spacing, int verletUpdateFrq);

    virtual void
    calculateForces(const int id_i, const int i);

    virtual void
    calculateStress(const int id_i, const int i,
                                 const int (&indexStress)[6]);
    virtual void
    updateState();

    void
    setForceScaling(double forceScaling);

    void
    setVerletRadius(double verletSpacing);

    virtual void
    initialize(double E, double nu, double delta, int dim, double h, double lc);

    virtual void
    applySurfaceCorrectionStep1(double strain);

    virtual void
    applySurfaceCorrectionStep2();
};
//------------------------------------------------------------------------------
}
#endif // CONTACTFORCE_H
