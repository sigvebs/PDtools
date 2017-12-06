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
void ContactForce::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    updateVerletList("contectForce", m_particles, m_grid, m_verletRadius);
    m_dim = dim;
}
//------------------------------------------------------------------------------
void ContactForce::applySurfaceCorrectionStep1(double strain)
{
    (void) strain;
}
//------------------------------------------------------------------------------
void ContactForce::applySurfaceCorrectionStep2()
{

}
//------------------------------------------------------------------------------
ContactForce::ContactForce(PD_Particles &particles, Grid & grid, double spacing, int verletUpdateFrq):
    Force(particles),
    m_grid(grid),
    m_steps(0),
    m_spacing(spacing),
    m_verletUpdateFrq(verletUpdateFrq)
{
    m_calulateStress = true;
    m_verletRadius = 2.0*spacing;
    m_forceScaling = 1.0;
    m_forceScaling = 15.0;

    m_forceScaling = 3.0;
    m_scaling = 0.95;
    m_scaling_dr0 = 0.9;

    m_scaling = 0.90;
    m_scaling_dr0 = 0.70;

    m_idToCol = &m_particles.idToCol();
    m_verletListId = particles.registerVerletList("contectForce");
    particles.setVerletUpdateFreq(10);
    m_indexMicromodulus = m_particles.registerParameter("micromodulus");
    m_indexRadius = m_particles.getParamId("radius");
    m_indexVolume = m_particles.getParamId("volume");
    m_iDr0 = m_particles.getPdParamId("dr0");
    m_indexConnected = m_particles.getPdParamId("connected");
    m_ghostParameters = {"volume", "micromodulus", "radius"};

    m_velocityScaling = 0.9999999999999999999999999995;
}
//------------------------------------------------------------------------------
void ContactForce::calculateForces(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);
    const double radius_i =  m_data(i, m_indexRadius);
    const vector<int> & verletList = m_particles.verletList(id_i);

    double dr_ij[m_dim];
    double dr0;

    bool contact = false;

    for(const int id_j:verletList) {
        const int j = (*m_idToCol).at(id_j);
        const double radius_j =  m_data(j, m_indexRadius);
//        double contactDistance = m_scaling*(radius_i + radius_j);
        double contactDistance = (radius_i + radius_j);
        double drSquared = 0;
        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            drSquared += dr_ij[d]*dr_ij[d];
        }

        if(drSquared < 1.8*contactDistance*contactDistance) {
            // Searching in the pd-connections
            bool hasBeenPdConnected = false;
            bool isPdConnected = false;
            const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

            for(const auto &con:PDconnections) {
                if(con.first == id_j) {
                    dr0 = con.second[m_iDr0];
                    if(con.second[m_indexConnected]) {
                        isPdConnected = true;
                    }
                    hasBeenPdConnected = true;
                    continue;
                }
            }
            if(isPdConnected)
                continue;

            if(hasBeenPdConnected) {
                contactDistance = min(contactDistance*m_scaling, dr0*m_scaling_dr0);
            }
            if(drSquared > contactDistance*contactDistance)
                continue;

            const double vol_j = m_data(j, m_indexVolume);
            const double c_j = m_data(j, m_indexMicromodulus);
            const double c_ij = 0.5*(c_i + c_j);

            double drLen = sqrt(drSquared);
            // To avoid roundoff errors
            if (drLen < 0)
                continue;

            contact = true;
            const double ds = (drLen - contactDistance)/contactDistance;
            const double fbond = m_forceScaling*c_ij*ds*vol_j/drLen;

            for(int d=0; d<m_dim; d++) {
                m_F(i, d) += dr_ij[d]*fbond;
                m_v *= m_velocityScaling;
            }
        }
    }

    if(contact) {
        for(int d=0; d<m_dim; d++) {
            m_v *= m_velocityScaling;
        }
    }
}
//------------------------------------------------------------------------------
void ContactForce::calculateStress(const int id_i, const int i,
                                   const int (&indexStress)[6])
{
    (void) id_i;
    (void) i;
    (void) indexStress;
    /*
    const double c_i = m_data(i, m_indexMicromodulus);
    const double radius_i =  m_data(i, m_indexRadius);

    double dr_ij[m_dim];

    const vector<int> & verletList = m_particles.verletList(id_i);

    for(const int id_j:verletList) {
        const int j = (*m_idToCol).at(id_j);
        const double radius_j =  m_data(i, m_indexRadius);
        const double contactDistance = m_scaling*(radius_i + radius_j);
        double drSquared = 0;
        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            drSquared += dr_ij[d]*dr_ij[d];
        }
        double drLen = sqrt(drSquared);

        if(drLen < contactDistance) {
            // Searching in the pd-connections
            bool isPdConnected = false;
            const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
            for(const auto &con:PDconnections) {
                if(con.first == id_j && con.second[m_indexConnected]) {
                    isPdConnected = true;
                    continue;
                }
            }
            if(isPdConnected)
                continue;

            const double vol_j = m_data(j, m_indexVolume);
            const double c_j = m_data(j, m_indexMicromodulus);
            const double c_ij = 0.5*(c_i + c_j);
            const double ds = (drLen - contactDistance)/contactDistance;
            const double bond_ij = m_forceScaling*c_ij*ds*vol_j/drLen;

            m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
            m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
            m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

            if(m_dim == 3) {
                m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
                m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
                m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
            }
        }
    }
    */
}
//------------------------------------------------------------------------------
void ContactForce::updateState()
{
    if (m_steps % m_verletUpdateFrq == 0) {
        updateVerletList("contectForce", m_particles, m_grid, m_verletRadius);
    }
    m_steps++;
}
//------------------------------------------------------------------------------
}

