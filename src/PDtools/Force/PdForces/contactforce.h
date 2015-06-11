#ifndef CONTACTFORCE_H
#define CONTACTFORCE_H

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
    arma::mat & m_F;
    arma::mat & m_r;
    arma::mat & m_data;
    double m_scaling;
    double m_verletRadius;
    double m_forceScaling;
    int m_indexVolume;
    int m_verletUpdateFrq = 20;
    int m_verletListId;
    int m_indexMicromodulus;
    int m_indexRadius;
public:
    ContactForce(PD_Particles &particles, Grid &grid, double spacing);
    ~ContactForce();
    virtual void calculateForces(const std::pair<int, int> & idCol);
    virtual void calculateStress(const std::pair<int, int> & idCol,
                                 const int (&indexStress)[6]);
    virtual void updateState();
    void setForceScaling(double forceScaling);
    void setVerletRadius(double verletSpacing);
};
//------------------------------------------------------------------------------
}
#endif // CONTACTFORCE_H
