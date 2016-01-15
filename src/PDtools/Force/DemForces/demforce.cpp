#include "demforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
DemForce::DemForce(PD_Particles &particles, double T, double C):
    Force(particles)
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexVolume = m_particles.getParamId("volume");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexStretch = m_particles.registerPdParameter("stretch");
    m_indexConnected = m_particles.getPdParamId("connected");
    m_indexRadius = m_particles.getParamId("radius");
    m_indexUnbreakable =  m_particles.registerParameter("unbreakable");

    m_ghostParameters = {"volume", "radius"};


    if(m_dim == 2)
    {
        m_indexFs[0] = m_particles.registerPdParameter("Fs_x");
        m_indexFs[1] = m_particles.registerPdParameter("Fs_y");
    }
    else if(m_dim == 3)
    {
        m_indexFs[0] = m_particles.registerPdParameter("Fs_x");
        m_indexFs[1] = m_particles.registerPdParameter("Fs_y");
        m_indexFs[2] = m_particles.registerPdParameter("Fs_z");
    }

    m_T = T;
    m_C = C;
}
//------------------------------------------------------------------------------
void DemForce::calculateForces(const int id_i, const int i)
{
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];
    const double radius_i = m_data(i, m_indexRadius);
    const double A_i = M_PI*pow(radius_i, 2);

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }
        const double dr = sqrt(dr2);
        const double dr0 = con.second[m_indexDr0];
        const double normalForce = m_Kn*(dr - dr0);
        const double bond_ij = normalForce/dr;

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*bond_ij;
        }

        const double radius_j = m_data(j, m_indexRadius);
        const double A_j = M_PI*pow(radius_j, 2);
        const double A = min(A_i, A_j);

        if(normalForce > m_T*A)
        {
            if(!m_data(j, m_indexUnbreakable))
                con.second[m_indexConnected] = 0;
        }
    }
}
//------------------------------------------------------------------------------
void DemForce::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];
    const double vol_i = m_data(i, m_indexVolume);

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }
        const double dr = sqrt(dr2);
        const double dr0 = con.second[m_indexDr0];
        const double normalForce = m_Kn*(dr - dr0);
        const double bond_ij = normalForce/dr/vol_i;

        m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
        m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
        }
    }
}
//------------------------------------------------------------------------------
void DemForce::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);

    double c, k;
    double A;
    if(dim == 3)
    {
        A = pow(lc, 3);
        nu = 1./4.;
        k = E/(3.*(1. - 2.*nu));
        c = 18.*k/(M_PI*pow(delta, 4))*pow(A,2);
    }
    else if(dim == 2)
    {
        A = pow(lc, 2)*h;
        nu = 1./3.;
        k = E/(2.*(1. - nu));
        c = 12*k/(h*M_PI*pow(delta, 3))*pow(A,2);
    }
    else if(dim == 1)
    {
        A = lc*pow(h, 2);
        nu = 1./4.;
        k = E;
        c = 2*E/(h*h*pow(delta, 2))*pow(A,2);
    }


    m_Kn = c;
    m_Ks = 0.35*c;
    cout << "c = " << m_Kn << endl;
}
//------------------------------------------------------------------------------
}
