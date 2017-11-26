#include "pd_lpsdampenedcontact_p.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_lpsDampenedContact_porosity::PD_lpsDampenedContact_porosity(PD_Particles &particles, double m, double b, double c, bool planeStress):
    PD_LPS_POROSITY(particles, m, b, planeStress),
    m_dampCoeff(c)
{
    particles.setNeedGhostVelocity(true);
}
//------------------------------------------------------------------------------
void PD_lpsDampenedContact_porosity::calculateForces(const int id, const int i)
{
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);
    const double a_i = m_data(i, m_iA);
    const double b_i = m_data(i, m_iA);

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
        const double a_j = m_data(j, m_iA);
        const double b_j = m_data(j, m_iB);

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
        double bond = (b_i*theta_i/m_i + b_j*theta_j/m_j)*dr0;
        bond += (a_i/m_i + a_j/m_j)*ds;
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
        m_data(i, m_iThetaNew) = m_dim/m_i*thetaNew;
}
//------------------------------------------------------------------------------
double PD_lpsDampenedContact_porosity::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double theta_i = this->computeDilation(id_i, i);
    double dr_ij[m_dim];
    const double a_i = m_data(i, m_iA);
    const double k_i = m_data(i, m_iK);

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
        const double a_j = m_data(j, m_iA);
        double dr2 = 0;
        double drdv = 0;

        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double dsdt = drdv/(dr);
        const double ds = dr - dr0 + m_dampCoeff*dsdt;
        const double extension_term =  (a_i + a_j)*w*(pow(ds - theta_i*dr0/m_dim, 2));

        W_i += (extension_term)*vol_j*volumeScaling;
    }
    W_i += k_i*(pow(theta_i, 2));

    return 0.5*W_i;
}
//------------------------------------------------------------------------------
}
