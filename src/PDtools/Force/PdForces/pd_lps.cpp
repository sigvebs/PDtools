#include "pd_lps.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_LPS::PD_LPS(PD_Particles &particles, bool planeStress):
    Force(particles),
    m_planeStress(planeStress)
{
    m_iMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_iTheta = m_particles.registerParameter("theta", 0);
    m_iThetaNew = m_particles.registerParameter("thetaNew", 0);
    m_iMass = m_particles.registerParameter("LPS_mass", 1);

    m_iVolume = m_particles.getParamId("volume");
    m_iDr0 = m_particles.getPdParamId("dr0");
    m_iVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_iStretch = m_particles.registerPdParameter("stretch");
    m_iConnected = m_particles.getPdParamId("connected");
    m_indexBrokenNow = m_particles.registerParameter("brokenNow", 0);

//    m_iForceScalingDilation = m_particles.registerPdParameter("forceScalingDilation", 1.);
//    m_iForceScalingBond = m_particles.registerPdParameter("forceScalingBond", 1.);
    m_ghostParameters.push_back("volume");
    m_ghostParameters.push_back("theta");
    m_ghostParameters.push_back("LPS_mass");

    m_initialGhostParameters = {"volume", "theta", "LPS_mass"};
    m_hasUpdateState = true;
}
//------------------------------------------------------------------------------
void PD_LPS::calculateForces(const int id, const int i)
{
    const double theta_i = m_data(i, m_iTheta);
//    const double theta_i = this->computeDilation(id, i);
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    double thetaNew = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
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
        const double w = 1./dr0;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        double bond = m_c*(theta_i/m_i + theta_j/m_j)*dr0;
        bond += m_alpha*(1./m_i + 1./m_j)*ds;
        bond *= w*vol_j*volumeScaling/dr;
        thetaNew += w*dr0*ds*vol_j*volumeScaling;

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*bond;
        }

        con.second[m_iStretch] = ds/dr0;
    }

    if(nConnections <= 3)
        m_data(i, m_iThetaNew) = 0;
    else
        m_data(i, m_iThetaNew) = m_t*thetaNew/m_i;
}
//------------------------------------------------------------------------------
double PD_LPS::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double theta_i = this->computeDilation(id_i, i);
    double dr_ij[m_dim];

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();

    double W_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double w = 1./dr0;
        const double volumeScaling = con.second[m_iVolumeScaling];
        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        const double extension_term =  m_alpha*w*(pow(ds - theta_i*dr0/m_dim, 2));

        W_i += (extension_term)*vol_j*volumeScaling;
    }
    W_i += m_k*(pow(theta_i, 2));

    return 0.5*W_i;
}
//------------------------------------------------------------------------------
double PD_LPS::computeDilation(const int id_i, const int i)
{
    const double m_i = m_data(i, m_iMass);
    double dr_ij[m_dim];

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();

    double theta_i = 0;

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        double dr2 = 0;
        const double w = 1./dr0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        theta_i += w*dr0*ds*vol_j*volumeScaling;
    }

    if(nConnections <= 3)
    {
        theta_i = 0;
    }
    return theta_i*m_t/m_i;
}
//------------------------------------------------------------------------------
void PD_LPS::calculatePotentialEnergy(const int id_i, const int i, int indexPotential)
{
    double vol_i = m_data(i, m_iVolume);
    m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i)*vol_i;
}
//------------------------------------------------------------------------------
void PD_LPS::updateState(int id, int i)
{
    (void) id;
    m_data(i, m_iTheta) = m_data(i, m_iThetaNew);
}
//------------------------------------------------------------------------------
double PD_LPS::calculateStableMass(const int id_a, const int a, double dt)
{
    dt *= 1.1;

    const double vol_a = m_data(a, m_iVolume);
    const double m_a = m_data(a, m_iMass);
    const arma::mat & R0 = m_particles.r0();

    double m[m_dim];
    double dr0[m_dim];
    for(int d=0; d<m_dim; d++)
    {
        m[d] = 0;
    }

    const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_a);

    double k[m_dim];

    for(int i=0; i<m_dim; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            k[d] = 0;
        }

        for(auto &con:PDconnections)
        {
            if(con.second[m_iConnected] <= 0.5)
                continue;

            const int id_b = con.first;
            const int b = m_idToCol.at(id_b);

            for(int d=0; d<m_dim; d++)
            {
                dr0[d] = R0(a, d) - R0(b, d);
            }

            const double dr0Len = con.second[m_iDr0];
            const double m_b = m_data(b, m_iMass);
            const double vol_b = m_data(b, m_iVolume);
            const double volumeScaling = con.second[m_iVolumeScaling];
            const double Va = vol_a*volumeScaling;
            const double Vb = vol_b*volumeScaling;
            const double w = 1./dr0Len;

            // Check this
            const double dr0Len2 = pow(dr0Len, 2);
//            double C = m_dim*m_c*(Vb/pow(m_a,2) + Va/pow(m_b,2)) + m_alpha*(1./m_a + 1./m_b);
            double C =  m_alpha*(1./m_a + 1./m_b);
//            double C = m_dim*m_c*(1./pow(m_a,2) + 1./pow(m_b,2)) + m_alpha*(1./m_a + 1./m_b);
            C *= w*Vb/dr0Len2;

            double sum = 0;

            for(int j=0; j<m_dim; j++)
            {
                sum += fabs(dr0[j]);
            }

            k[i] += fabs(dr0[i])*C*sum;
        }
        m[i] = k[i];
    }

    double stiffness = 0;

    for(int d=0;d<m_dim; d++)
    {
        if(m[d]>stiffness)
        {
            stiffness = m[d];
        }
    }

    return 4.*0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void PD_LPS::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    m_delta = delta;
    m_dim = dim;
    m_nu = nu;
    m_mu = 0.5*E/(1 + nu);
    m_k = E/(3.*(1. - 2.*nu));

    double k  = E/(2.*(1. - nu));
    double c = 2.*(k - 2.*m_mu);

    if(dim == 3)
    {
        m_t = 3.;
        m_c = (3.*m_k - 5.*m_mu);
        m_alpha = 15.*m_mu;
    }
    else if(dim == 2)
    {
        m_alpha = 8.*m_mu;
        if(m_planeStress)
        {
            m_t = 2.*(2.*m_nu - 1.)/(m_nu - 1.);
            double k_ = m_k + m_mu/9.*pow((m_nu + 1.)/(2*m_nu - 1), 2);
            m_c = m_t*k_ - 8./3.*m_mu*(2. - m_t/3.);
        }
        else // Plane strain
        {
            m_t = 2.;
            m_c = 2.*(m_k - 15./9.*m_mu);
        }
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported." << endl;
        cerr << "use 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    calculateWeightedVolume();
/*
    if(dim == 3)
    {
        m_k = E/(3.*(1. - 2.*nu));
        m_c = (3.*m_k - 5.*m_mu);
        m_alpha = 15.*m_mu;
    }
    else if(dim == 2)
    {
        m_k  = E/(2.*(1. - nu));
        m_c = 2.*(m_k - 2.*m_mu);
        m_alpha = 8.*m_mu;
    }
    */
}
//------------------------------------------------------------------------------
void PD_LPS::calculateWeightedVolume()
{
    const int nParticles = m_particles.nParticles();
    bool analytical = false;

    // Calculating the one-body forces
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = m_colToId.at(i);

        const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
        const int nConnections = PDconnections.size();
        double m = 0;
        for(int l_j=0; l_j<nConnections; l_j++)
        {
            auto &con = PDconnections[l_j];
            if(con.second[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_idToCol.at(id_j);
            const double volumeScaling = con.second[m_iVolumeScaling];
            const double vol_j = m_data(j, m_iVolume);
            const double dr0 = con.second[m_iDr0];
            const double w = 1./dr0;

            m += w*dr0*dr0*vol_j*volumeScaling;
        }
        if(analytical)
        {
            if(m_dim == 3)
            {
                m = 4.*M_PI/5.*pow(m_delta, 5);
            }else
            {
                m = m_h*M_PI/2.*pow(m_delta, 4);
            }

        }

        m_data(i, m_iMass) = m;

        // Setting the micromodulus
        m_data(i, m_iMicromodulus) = m_alpha/m;
    }
}
//------------------------------------------------------------------------------
void PD_LPS::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{
    /*
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    double beta = 0;
    if(m_dim == 3)
        beta = 15*m_mu;
    else
        beta = 8*m_mu;

    const double alpha_i = beta/m_i;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    for(int l_j=0; l_j<nConnections; l_j++)
    {
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
        const double alpha_j = beta/m_j;
        const double gb_ij = con.second[m_iForceScalingBond];
        const double gd_ij = con.second[m_iForceScalingDilation];

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;

        double bond_ij = gd_ij*(3*m_k - 5*m_mu)*(theta_i/m_i + theta_j/m_j)*dr0;
        bond_ij += gb_ij*(alpha_i + alpha_j)*ds;
        bond_ij *= vol_j*volumeScaling/dr;

        m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
        m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
        }
    }
    */
}
//------------------------------------------------------------------------------
}
