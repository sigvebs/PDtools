#include "contactforce.h"

#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
void ContactForce::setForceScaling(double forceScaling)
{
    m_forceScaling = forceScaling;
}
//------------------------------------------------------------------------------
void ContactForce::setVerletRadius(double verletSpacing)
{
    m_verletRadius = verletSpacing;
}
//------------------------------------------------------------------------------
ContactForce::ContactForce(PD_Particles &particles, Grid & grid, double spacing):
    Force(particles),
    m_grid(grid),
    m_steps(0),
    m_F(m_particles.F()),
    m_r(m_particles.r()),
    m_data(m_particles.data())
{
    m_verletRadius = 1.5*spacing;
    m_forceScaling = 15.0;
    m_scaling = 0.48;

    m_verletListId = particles.registerVerletList("contectForce");
    particles.setVerletUpdateFreq(30);
    m_indexMicromodulus = m_particles.getParamId("micromodulus");
    m_indexRadius = m_particles.getParamId("radius");
    m_indexVolume = m_particles.getParamId("volume");
}
//------------------------------------------------------------------------------
ContactForce::~ContactForce()
{
}
//------------------------------------------------------------------------------
void ContactForce::calculateForces(const std::pair<int, int> &idCol)
{
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    double radius_i =  m_data(col_i, m_indexRadius);
    double x_i = m_r(X, col_i);
    double y_i = m_r(Y, col_i);
    double z_i = m_r(Z, col_i);
    double dr_ij[3];

    const vector<int> & verletList = m_particles.verletList(pId);

    for(int col_j:verletList)
    {
        double radius_j =  m_data(col_i, m_indexRadius);
        double contactDistance = m_scaling*(radius_i + radius_j);

        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);

        if(drLen < contactDistance)
        {
            double vol_j = m_data(col_j, m_indexVolume);
            double c_j = m_data(col_j, m_indexMicromodulus);
            double c_ij = 0.5*(c_i + c_j);
            double ds = drLen - contactDistance;
            double fbond = m_forceScaling*c_ij*ds*vol_j/drLen;

            m_F(X, col_i) += dr_ij[X]*fbond;
            m_F(Y, col_i) += dr_ij[Y]*fbond;
            m_F(Z, col_i) += dr_ij[Z]*fbond;
        }
    }
}
//------------------------------------------------------------------------------
void ContactForce::calculateStress(const std::pair<int, int> &idCol,
                                   const int (&indexStress)[6])
{
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    double x_i = m_r(X, col_i);
    double y_i = m_r(Y, col_i);
    double z_i = m_r(Z, col_i);
    double dr_ij[3];
    double f[3];
    double radius_i =  m_data(col_i, m_indexRadius);

    const vector<int> & verletList = m_particles.verletList(pId);

    for(int col_j:verletList)
    {
        double radius_j =  m_data(col_i, m_indexRadius);
        double contactDistance = m_scaling*(radius_i + radius_j);
        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);

        if(drLen < contactDistance)
        {
            double vol_j = m_data(col_j, m_indexVolume);
            double c_j = m_data(col_j, m_indexMicromodulus);
            double c_ij = 0.5*(c_i + c_j);
            double ds = drLen - contactDistance;
            double fbond = m_forceScaling*c_ij*ds*vol_j/drLen;

            f[X] = dr_ij[X]*fbond;
            f[Y] = dr_ij[Y]*fbond;
            f[Z] = dr_ij[Z]*fbond;

            m_data(col_i, indexStress[0]) += 0.5*f[X]*dr_ij[X];
            m_data(col_i, indexStress[1]) += 0.5*f[Y]*dr_ij[Y];
            m_data(col_i, indexStress[2]) += 0.5*f[Z]*dr_ij[Z];
            m_data(col_i, indexStress[3]) += 0.5*f[X]*dr_ij[Y];
            m_data(col_i, indexStress[4]) += 0.5*f[X]*dr_ij[Z];
            m_data(col_i, indexStress[5]) += 0.5*f[Z]*dr_ij[Z];
        }
    }
}
//------------------------------------------------------------------------------
void ContactForce::updateState()
{
    if (m_steps % m_verletUpdateFrq == 0)
    {
        updateVerletList("contectForce", m_particles, m_grid, m_verletRadius);
    }
    m_steps++;
}
//------------------------------------------------------------------------------
}

