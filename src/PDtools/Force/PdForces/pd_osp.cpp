#include "pd_osp.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_OSP::PD_OSP(PD_Particles &particles):
    Force(particles)
{
    m_indexA = m_particles.registerParameter("a", 1);
//    m_indexB = m_particles.registerParameter("b", 1);
    m_indexB = m_particles.registerParameter("micromodulus", 1);
    m_indexD = m_particles.registerParameter("d", 1);
    m_indexTheta = m_particles.registerParameter("theta", 0);
    m_indexThetaNew = m_particles.registerParameter("thetaNew", 0);

    m_indexVolume = m_particles.getParamId("volume");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexStretch = m_particles.registerPdParameter("stretch");
    m_indexConnected = m_particles.getPdParamId("connected");

    m_indexForceScalingDilation = m_particles.registerPdParameter("forceScalingDilation", 1.);
    m_indexForceScalingBond = m_particles.registerPdParameter("forceScalingBond", 1.);

    m_hasUpdateState = true;
}
//------------------------------------------------------------------------------
PD_OSP::~PD_OSP()
{

}
//------------------------------------------------------------------------------
void PD_OSP::calculateForces(const int id, const int i)
{
    const double a_i = m_data(i, m_indexA);
    const double b_i = m_data(i, m_indexB);
    const double d_i = m_data(i, m_indexD);
    const double theta_i = m_data(i, m_indexTheta);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    double f_i[3];
    for(int d=0; d<m_dim; d++)
    {
        f_i[d] = 0;
    }

    const int nConnections = PDconnections.size();
    double dr_ij[3];

    double thetaNew = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double a_j = m_data(j, m_indexA);
        const double b_j = m_data(j, m_indexB);
        const double d_j = m_data(j, m_indexD);
        const double theta_j = m_data(j, m_indexTheta);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double a_ij = 0.5*(a_i + a_j);
        const double b_ij = 0.5*(b_i + b_j);
        const double d_ij = 0.5*(d_i + d_j);
        const double Gd_ij = con.second[m_indexForceScalingDilation];
        const double Gb_ij = con.second[m_indexForceScalingBond];

        double dr2 = 0;
        double A_ij = 0; // The lambda-factor

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            A_ij += dr_ij[d]*(m_r0(j, d) - m_r0(i, d));
        }

        const double dr = sqrt(dr2);
        A_ij /= (dr0*dr);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        const double fbond = (a_ij*d_ij*Gd_ij*A_ij/dr0*(theta_i + theta_j) + b_ij*Gb_ij*s)
                *vol_j*volumeScaling/dr;
        for(int d=0; d<m_dim; d++)
        {
            f_i[d] += dr_ij[d]*fbond;
        }

        thetaNew += d_ij*s*A_ij*vol_j*volumeScaling;
        con.second[m_indexStretch] = s;
    }

    for(int d=0; d<m_dim; d++)
    {
        m_F(i, d) += m_delta*f_i[d];
    }

    m_data(i, m_indexThetaNew) = m_delta*thetaNew;
}
//------------------------------------------------------------------------------
double PD_OSP::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double a_i = m_data(i, m_indexA);
    const double b_i = m_data(i, m_indexB);
    const double d_i = m_data(i, m_indexD);
//    const double theta_i = m_data(i, m_indexTheta);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections.size();
    double dr_ij[3];

    double theta_sum = 0;
    double W_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double b_j = m_data(j, m_indexB);
        const double d_j = m_data(j, m_indexD);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double b_ij = 0.5*(b_i + b_j);
        const double d_ij = 0.5*(d_i + d_j);
        const double Gd_ij = con.second[m_indexForceScalingDilation];
        const double Gb_ij = con.second[m_indexForceScalingBond];

        double dr2 = 0;
        double A_ij = 0; // The lambda-factor

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            A_ij += dr_ij[d]*(m_r0(j, d) - m_r0(i, d));
        }

        const double dr = sqrt(dr2);
        A_ij /= (dr0*dr);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        theta_sum += d_ij*Gd_ij*A_ij*s*vol_j*volumeScaling;
        W_i += b_ij*Gb_ij*s*s*dr*vol_j*volumeScaling;
    }

    return m_delta*(a_i*theta_sum*theta_sum + W_i);
}
//------------------------------------------------------------------------------
void PD_OSP::calculatePotentialEnergy(const int id_i, const int i, int indexPotential)
{
    double vol_i = m_data(i, m_indexVolume);
    m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i)*vol_i;
}
//------------------------------------------------------------------------------
void PD_OSP::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{
    const double a_i = m_data(i, m_indexA);
    const double b_i = m_data(i, m_indexB);
    const double d_i = m_data(i, m_indexD);
    const double theta_i = m_data(i, m_indexTheta);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections.size();
    double dr_ij[3];
    double f[3];

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double a_j = m_data(j, m_indexA);
        const double b_j = m_data(j, m_indexB);
        const double d_j = m_data(j, m_indexD);
        const double theta_j = m_data(j, m_indexTheta);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double a_ij = 0.5*(a_i + a_j);
        const double b_ij = 0.5*(b_i + b_j);
        const double d_ij = 0.5*(d_i + d_j);
        const double Gd_ij = con.second[m_indexForceScalingDilation];
        const double Gb_ij = con.second[m_indexForceScalingBond];

        double dr2 = 0;
        double A_ij = 0; // The lambda-factor

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            A_ij += dr_ij[d]*(m_r0(j, d) - m_r0(i, d));
        }

        const double dr = sqrt(dr2);
        A_ij /= (dr0*dr);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        const double fbond = m_delta*(a_ij*d_ij*Gd_ij*A_ij/dr0*(theta_i + theta_j) + b_ij*Gb_ij*s)
                *vol_j*volumeScaling/dr;

        for(int d=0; d<m_dim; d++)
        {
            f[d] = dr_ij[d]*fbond;
        }

        m_data(i, indexStress[0]) += 0.5*f[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*f[Y]*dr_ij[Y];
        m_data(i, indexStress[2]) += 0.5*f[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[3]) += 0.5*f[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*f[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*f[Y]*dr_ij[Z];
        }
    }
}
//------------------------------------------------------------------------------
void PD_OSP::updateState(int id, int i)
{
    (void) id;
    m_data(i, m_indexTheta) = m_data(i, m_indexThetaNew);
}
//------------------------------------------------------------------------------
void PD_OSP::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    double a, b, d;
    double k; // Bulk modulus
    double mu; // Shear modulus

    m_delta = delta;
    m_dim = dim;
    mu = 0.5*E/(1 + nu);

    if(dim == 3)
    {
        k = E/(3.*(1. - 2.*nu));
        a = 0.5*(k - 5.*mu/3.);
        b = 15.*mu/(2*M_PI*pow(delta, 5));
        d = 9./(4.*M_PI*pow(delta, 4));
    }
    else if(dim == 2)
    {
        k = E/(2.*(1. - nu));
        a = 0.5*(k - 2*mu);
        b = 6.*mu/(h*M_PI*pow(delta, 4));
        d = 2./(h*M_PI*pow(delta, 3));
    }
    else if(dim == 1)
    {
        double A = h*h;
        a = 0;
        b = E/(2.*A*pow(delta, 3));
        d = 1.0/(2.*pow(delta, 2)*A);
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported." << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    m_particles.setParameter("a", a);
//    m_particles.setParameter("b", b);
    m_particles.setParameter("micromodulus", b);
    m_particles.setParameter("d", d);

    if(m_numericalInitialization){
        double strain = 0.001;
        applySurfaceCorrectionStep1(mu, nu, dim, strain);
    }
}
//------------------------------------------------------------------------------
void PD_OSP::applySurfaceCorrectionStep1(double mu, double nu, int dim, double strain)
{
    (void) nu;
    arma::vec3 strainFactor;
    arma::mat gd = arma::zeros(m_particles.nParticles(), dim); // Dilation correction
    arma::mat gb = arma::zeros(m_particles.nParticles(), dim); // Bond correction

    cerr << "TODO: fix the PD_OSP::applySurfaceCorrectionStep1" << endl;
    exit(1);
    /*
    //--------------------------------------------------------------------------
    // Apllying correction to the dilation term
    //--------------------------------------------------------------------------
    // Stretching all particle in the x, y and z-direction
    strainFactor(0) = strain;
    strainFactor(1) = 0;
    strainFactor(2) = 0;

    for(int a=0; a<dim; a++)
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
            int col_i = idCol.second;

            for(int d=0; d<dim; d++)
            {
                m_r(col_i, d) = (1 + strainFactor(d))*m_r(col_i, d);
            }
        }

        double dilation_term = strain;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double theta_i = calculateDilationTerm(idCol);
            double factor =  dilation_term/theta_i;
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

            for(int d=0; d<dim; d++)
            {
                m_r(col_i, d) = m_r(col_i, d)/(1 + strainFactor(d));
            }
        }
    }

    //--------------------------------------------------------------------------
    // Applying correction to the bond term
    //--------------------------------------------------------------------------
    // Performing a simple shear of all particle in the x, y and z-direction
    arma::ivec3 axis;
    strainFactor(0) = 0.5*strain;
    strainFactor(1) = 0.*strain;
    strainFactor(2) = 0;
    axis(0) = 1;
    axis(1) = 0;
    axis(2) = 0;

    for(int a=0; a<dim; a++)
    {
        if(a == 1)
        {
            strainFactor.swap_rows(1,2);
            strainFactor.swap_rows(0,1);
            axis(0) = 0;
            axis(1) = 2;
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

            for(int d=0; d<dim; d++)
            {
                double shear = strainFactor(d)*m_r(axis(d), col_i);
                m_r(col_i, d) = m_r(col_i, d) + shear;
//                m_r(col_i, d) = m_r(col_i, d) + strainFactor(d)*m_r(axis(d), col_i);
//                m_r(col_i, d) = (1 + strainFactor(d))*m_r(col_i, d);
            }
        }

        double dilation_term = 0.5*mu*strain*strain;
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;
            double bond_i = calculateBondPotential(idCol);
            double factor =  dilation_term/bond_i;
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

            for(int d=0; d<dim; d++)
            {
                m_r(col_i, d) = m_r(col_i, d) - strainFactor(d)*m_r(axis(d), col_i);
//                m_r(col_i, d) = m_r(col_i, d)/(1 + strainFactor(d));
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
            int col_j = m_idToCol.at(id_j);

            double dr0Len = con.second[m_indexDr0];
            arma::vec3 n = (m_r.col(col_i) - m_r.col(col_j))/dr0Len;

            arma::vec3 gd_mean;
            arma::vec3 gb_mean;
            double Gd = 0;
            double Gb = 0;
            for(int d=0; d<dim; d++)
            {
                gd_mean(d) = 0.5*(gd(col_i, d) + gd(col_j, d));
                gb_mean(d) = 0.5*(gb(col_i, d) + gb(col_j, d));
                Gd += pow(n(d)/gd_mean(d), 2);
                Gb += pow(n(d)/gb_mean(d), 2);
            }

            Gd = pow(Gd, -0.5);
            Gb = pow(Gb, -0.5);
            con.second[m_indexForceScalingDilation] *= Gd;
            con.second[m_indexForceScalingBond] *= Gb;
        }
    }
    */
}
//------------------------------------------------------------------------------
double PD_OSP::calculateDilationTerm(const int id_i, const int i)
{
    const double d_i = m_data(i, m_indexD);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const int nConnections = PDconnections.size();
    double dr_ij[3];

    double theta_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double d_j = m_data(j, m_indexD);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double d_ij = 0.5*(d_i + d_j);

        double dr2 = 0;
        double A_ij = 0; // The lambda-factor

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            A_ij += dr_ij[d]*(m_r0(d, j) - m_r0(i, d));
        }

        const double dr = sqrt(dr2);
        A_ij /= (dr0*dr);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        theta_i += d_ij*s*A_ij*vol_j*volumeScaling;
    }

    return m_delta*theta_i;
}
//------------------------------------------------------------------------------
double PD_OSP::calculateBondPotential(const int id_i, const int i)
{
    const double b_i = m_data(i, m_indexB);
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections.size();
    double dr_ij[3];

    double bond_i = 0;
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double b_j = m_data(j, m_indexD);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double b_ij = 0.5*(b_i + b_j);

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        bond_i += b_ij*dr0*s*s*vol_j*volumeScaling;
    }

    return m_delta*bond_i;
}
//------------------------------------------------------------------------------
}
