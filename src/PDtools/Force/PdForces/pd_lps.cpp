#include "pd_lps.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_LPS::PD_LPS(PD_Particles &particles):
    Force(particles)
{
    m_iTheta = m_particles.registerParameter("theta", 0);
    m_iThetaNew = m_particles.registerParameter("thetaNew", 0);
    m_iMass = m_particles.registerParameter("LMP_mass", 1);

    m_iVolume = m_particles.getParamId("volume");
    m_iDr0 = m_particles.getPdParamId("dr0");
    m_iVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_iStretch = m_particles.registerPdParameter("stretch");
    m_iConnected = m_particles.getPdParamId("connected");

    m_iForceScalingDilation = m_particles.registerPdParameter("forceScalingDilation", 1.);
    m_iForceScalingBond = m_particles.registerPdParameter("forceScalingBond", 1.);
    m_iCompute = m_particles.getPdParamId("compute");
}
//------------------------------------------------------------------------------
PD_LPS::~PD_LPS()
{

}
//------------------------------------------------------------------------------
void PD_LPS::calculateForces(const std::pair<int, int> &idCol)
{
    const int pId = idCol.first;
    const int i = idCol.second;
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    double alpha;
    if(m_dim == 3)
        alpha = 15*m_mu;
    else
        alpha = 8*m_mu;

    const double c = (3*m_k - 5*m_mu);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    double thetaNew = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];

        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double m_j = m_data(j, m_iMass);
        const double theta_j = m_data(j, m_iTheta);
        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
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
        double bond = gd_ij*c*(theta_i/m_i + theta_j/m_j)*dr0;
        bond += gb_ij*alpha*(1./m_i + 1./m_j)*ds;
        bond *= vol_j*volumeScaling/dr;
        thetaNew += dr0*ds*vol_j*volumeScaling;

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*bond;
        }

        con.second[m_iStretch] = ds/dr0;
    }

    m_data(i, m_iThetaNew) = m_dim/m_i*thetaNew;
}
//------------------------------------------------------------------------------
double PD_LPS::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    const int pId = idCol.first;
    const int i = idCol.second;
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    double alpha = 0;
    if(m_dim == 3)
        alpha = 15*m_mu/m_i;
    else
        alpha = 8*m_mu/m_i;

    double theta_i = this->computeDilation(idCol);

    double W_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        double ds = dr - dr0;
        const double extension_term =  alpha*(pow(ds - theta_i*dr0/m_dim, 2));

        W_i += (extension_term)*vol_j*volumeScaling;
    }
    W_i += m_k*(pow(theta_i, 2));

    return 0.5* W_i;
}
//------------------------------------------------------------------------------
double PD_LPS::computeDilation(const std::pair<int, int> &idCol)
{
    const int pId = idCol.first;
    const int i = idCol.second;
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    const int nConnections = PDconnections.size();

    double dr_ij[m_dim];
    double theta_i = 0;

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        theta_i += dr0*ds*vol_j*volumeScaling;
    }

    return theta_i*m_dim/m_i;
}
//------------------------------------------------------------------------------
void PD_LPS::calculatePotentialEnergy(const std::pair<int, int> &idCol, int indexPotential)
{
    int col_i = idCol.second;
    double vol_i = m_data(col_i, m_iVolume);
    m_data(col_i, indexPotential) += calculatePotentialEnergyDensity(idCol)*vol_i;
}
//------------------------------------------------------------------------------
void PD_LPS::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    const int pId = idCol.first;
    const int i = idCol.second;
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    double beta = 0;
    if(m_dim == 3)
        beta = 15*m_mu;
    else
        beta = 8*m_mu;

    const double alpha_i = beta/m_i;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

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
void PD_LPS::updateState(const std::pair<int, int> &idCol)
{
    const int col_i = idCol.second;
    m_data(col_i, m_iTheta) = m_data(col_i, m_iThetaNew);
}
//------------------------------------------------------------------------------
double PD_LPS::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    const int pId = idCol.first;
    const int a = idCol.second;
    const double vol_a = m_data(a, m_iVolume);
    const double m_a = m_data(a, m_iMass);
    const arma::mat & matR0 = m_particles.r0();

    double beta = 0;
    if(m_dim == 3)
        beta = 15*m_mu;
    else
        beta = 8*m_mu;

    const double alpha_a = beta/m_a;

    double m[m_dim];
    double dr0[m_dim];
    for(int d=0; d<m_dim; d++)
    {
        m[d] = 0;
    }

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

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
            const int b = m_pIds[id_b];

            for(int d=0; d<m_dim; d++)
            {
                dr0[d] = matR0(a, d) - matR0(b, d);
            }

            const double m_b = m_data(b, m_iMass);
            const double dr0Len = con.second[m_iDr0];
            const double vol_b = m_data(b, m_iVolume);
            const double volumeScaling = con.second[m_iVolumeScaling];
            const double Va = vol_a *volumeScaling;
            const double Vb = vol_b *volumeScaling;
            const double alpha_b = beta/m_b;
            const double gb_ij = con.second[m_iForceScalingBond];
            const double gd_ij = con.second[m_iForceScalingDilation];

            const double dr0Len2 = pow(dr0Len, 2);
            double coeff = gd_ij*3*dr0Len*(3*m_k - 5*m_mu)*(Vb/pow(m_a,2) + Va/pow(m_b,2));
            coeff += gb_ij*(alpha_a + alpha_b);
            coeff *= Vb/dr0Len2;
            double sum_xyz = 0;

            for(int j=0; j<m_dim; j++)
            {
                sum_xyz += fabs(dr0[j]);
            }

            k[i] += fabs(dr0[i])*coeff*sum_xyz;
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
    double stableMass = 4*0.25*pow(dt, 2)*stiffness;
    return stableMass;
}
//------------------------------------------------------------------------------
void PD_LPS::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    m_delta = delta;
    m_dim = dim;
    m_mu = 0.5*E/(1 + nu);
    m_nu = nu;

    if(dim == 3)
    {
        m_k = E/(3.*(1. - 2.*nu));
    }
    else if(dim == 2)
    {
        m_k  = E/(2.*(1. - nu));
    }
    else if(dim == 1)
    {
        m_k = E;
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported." << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }
    calculateMass();
}
//------------------------------------------------------------------------------
void PD_LPS::calculateMass()
{
    bool analytical = false;

    // Calculateing the one-body forces
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        pair<int, int> id_col(i, i);

        const int pId = id_col.first;
        const int col_i = id_col.second;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        const int nConnections = PDconnections.size();
        double m = 0;
        for(int l_j=0; l_j<nConnections; l_j++)
        {
            auto &con = PDconnections[l_j];
            if(con.second[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_pIds[id_j];
            const double volumeScaling = con.second[m_iVolumeScaling];
            const double vol_j = m_data(j, m_iVolume);
            const double dr0 = con.second[m_iDr0];

            m += dr0*dr0*vol_j*volumeScaling;
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
        m_data(col_i, m_iMass) = m;
    }
}
//------------------------------------------------------------------------------
void PD_LPS::applySurfaceCorrection(double strain)
{
    arma::vec3 strainFactor;
    arma::mat gd = arma::zeros(m_particles.nParticles(), m_dim); // Dilation correction
    arma::mat gb = arma::zeros(m_particles.nParticles(), m_dim); // Bond correction

    //--------------------------------------------------------------------------
    // Apllying correction to the dilation term
    //--------------------------------------------------------------------------
    // Stretching all particle in the x, y and z-direction
    strainFactor(0) = strain;
    strainFactor(1) = 0;
    strainFactor(2) = 0;

    for(int a=0; a<m_dim; a++)
    {
        if(a == 1)
            strainFactor.swap_rows(0,1);
        else if(a == 2)
            strainFactor.swap_rows(1,2);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Applying uniaxial stretch
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                m_r(col_i, d) = (1 + strainFactor(d))*m_r(col_i, d);
            }
        }

        double W_d = strain;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            const double theta_i = computeDilation(idCol);
            const double factor =  W_d/theta_i;
            gd(col_i, a) = factor;
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Resetting the positions
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
//                m_r(col_i, d) = m_r(col_i, d)/(1 + strainFactor(d));
                m_r(col_i, d) = m_r0(col_i, d);
            }
        }
    }

    //--------------------------------------------------------------------------
    // Applying correction to the bond term
    //--------------------------------------------------------------------------
    // Performing a simple shear of all particle in the x, y and z-direction
    arma::ivec3 axis;
    strainFactor(0) = strain;
    strainFactor(1) = 0.*strain;
    strainFactor(2) = 0;
    axis(0) = 1;
    axis(1) = 0;
    axis(2) = 0;

    for(int a=0; a<m_dim; a++)
    {
        if(a == 1)
        {
            strainFactor.swap_rows(1,2);
            strainFactor.swap_rows(0,1);
            axis(0) = 2;
            axis(1) = 0;
            axis(2) = 1;
        }
        else if(a == 2)
        {
            strainFactor.swap_rows(2,0);
            strainFactor.swap_rows(1,2);
            axis(0) = 2;
            axis(1) = 0;
            axis(2) = 0;
        }
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Applying uniaxial stretch
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                double shear = strainFactor(d)*m_r(axis(d), col_i);
                m_r(col_i, d) = m_r(col_i, d) + shear;
            }
        }

        double W_s = 0.5*m_mu*strain*strain;
//#ifdef USE_OPENMP
        //# pragma omp parallel for
        //#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            const double bond_i = calculatePotentialEnergyDensity(idCol);
            const double factor =  W_s/bond_i;
            gb(col_i, a) = factor;
        }

        //#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Resetting the positions
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
//                m_r(col_i, d) = m_r(col_i, d) - strainFactor(d)*m_r(axis(d), col_i);
                m_r(col_i, d) = m_r0(col_i, d);
            }
        }
    }

    //--------------------------------------------------------------------------
    // Calculating the scaling
    //--------------------------------------------------------------------------
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        int pId = idCol.first;
        int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];

            double dr0Len = con.second[m_iDr0];
            arma::vec3 n = (m_r.row(col_i).t() - m_r.row(col_j).t())/dr0Len;

            arma::vec3 gd_mean;
            arma::vec3 gb_mean;
            double Gd = 0;
            double Gb = 0;
            for(int d=0; d<m_dim; d++)
            {
                const double gd_i = gd(col_i, d);
                const double gd_j = gd(col_j, d);
                const double gb_i = gb(col_i, d);
                const double gb_j = gb(col_j, d);
                gd_mean(d) = 0.5*(gd_i + gd_j);
                gb_mean(d) = 0.5*(gb_i + gb_j);
                Gb += pow(n(d)/gb_mean(d), 2);
                Gd += pow(n(d)/gd_mean(d), 2);
            }

            Gd = pow(Gd, -0.5);
            Gb = pow(Gb, -0.5);
            con.second[m_iForceScalingDilation] *= Gd;
            con.second[m_iForceScalingBond] *= Gb;
        }
    }
}
//------------------------------------------------------------------------------
}
