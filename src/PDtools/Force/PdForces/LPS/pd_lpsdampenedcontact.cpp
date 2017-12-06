#include "pd_lpsdampenedcontact.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_lpsDampenedContact::PD_lpsDampenedContact(PD_Particles &particles, double c, bool planeStress):
    PD_LPS(particles, planeStress),
    m_dampCoeff(c)
{
    particles.setNeedGhostVelocity(true);
}
//------------------------------------------------------------------------------
void PD_lpsDampenedContact::calculateForces(const int id, const int i)
{
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    double thetaNew = 0;
    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con = PDconnections[l_j];

        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double m_j = m_data(j, m_iMass);
        const double theta_j = m_data(j, m_iTheta);
        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        const double w = weightFunction(dr0);

        double dr2 = 0;
        double drdv = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double dsdt = drdv/(dr);
        const double ds = dr - dr0 + m_dampCoeff*dsdt;
        double bond = m_c*(theta_i*m_i + theta_j*m_j)*dr0;
        bond += m_alpha*(m_i + m_j)*ds;
        bond *= w*vol_j*volumeScaling/dr;
        thetaNew += w*dr0*ds*vol_j*volumeScaling;

        for(int d=0; d<m_dim; d++) {
            m_F(i, d) += dr_ij[d]*bond;
        }

        con.second[m_iStretch] = ds/dr0;
    }

    if(nConnections <= 3)
        m_data(i, m_iThetaNew) = 0;
    else
        m_data(i, m_iThetaNew) = m_dim*m_i*thetaNew;
}
//------------------------------------------------------------------------------
double PD_lpsDampenedContact::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double theta_i = this->computeDilation(id_i, i);
    double dr_ij[m_dim];

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();

    double W_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double w = weightFunction(dr0);
        const double volumeScaling = con.second[m_iVolumeScaling];
        double dr2 = 0;
        double drdv = 0;

        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double dsdt = drdv/(dr);
        const double ds = dr - dr0 + m_c*dsdt;
        const double extension_term =  m_alpha*w*(pow(ds - theta_i*dr0/m_dim, 2));

        W_i += (extension_term)*vol_j*volumeScaling;
    }
    W_i += m_k*(pow(theta_i, 2));

    return 0.5*W_i;
}
//------------------------------------------------------------------------------
}
