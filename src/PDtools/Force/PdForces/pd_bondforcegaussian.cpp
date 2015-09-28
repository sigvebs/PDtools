#include "pd_bondforcegaussian.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondforceGaussian::PD_bondforceGaussian(PD_Particles &particles,
                                           std::string weightType):
    Force(particles),
    m_r(m_particles.r()),
    m_F(m_particles.F()),
    m_data(m_particles.data()),
    m_pIds(m_particles.pIds()),
    m_weightType(weightType)
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexWeightFunction = m_particles.registerPdParameter("weightFunction", 1);
    if(m_particles.hasParameter("micromodulus"))
    {
        m_indexMicromodulus = m_particles.getParamId("micromodulus");
    }
    else
    {
        cerr << "ERROR: Particles data does not contain either"
             << " 'micromodulus'. This is needed for the Bondforce." << endl
             << "Errorcode: " << MicrmodulusNotSet << endl;
        throw MicrmodulusNotSet;
    }
    m_indexVolume = m_particles.getParamId("volume");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexForceScaling = m_particles.registerPdParameter("forceScalingBond", 1);
    m_indexStretch = m_particles.registerPdParameter("stretch");
    m_indexConnected = m_particles.getPdParamId("connected");
}
//------------------------------------------------------------------------------
PD_bondforceGaussian::~PD_bondforceGaussian()
{

}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculateForces(const std::pair<int, int> &idCol)
{
    // PD_bond
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double dr_ij[m_dim];

    for(auto &con:PDconnections)
    {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int col_j = m_pIds[id_j];

        const double c_j = m_data(col_j, m_indexMicromodulus);
        const double vol_j = m_data(col_j, m_indexVolume);
        const double dr0         = con.second[m_indexDr0];
        const double volumeScaling   = con.second[m_indexVolumeScaling];
        const double g_ij = con.second[m_indexForceScaling];
        const double w_ij = con.second[m_indexWeightFunction];
        const double c_ij = 0.5*(c_i + c_j)*w_ij*g_ij;

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(d, col_j) - m_r(d, i);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        const double s = ds/dr0;
        const double fbond_ij = c_ij*s*vol_j*volumeScaling/dr;
#ifdef USE_N3L
        const double fbond_ji = -c_ij*s*vol_i*volumeScaling/dr;
#endif

        for(int d=0; d<m_dim; d++)
        {
            m_F(d, i) += dr_ij[d]*fbond_ij;
#ifdef USE_N3L
            m_F(d, j) += dr_ij[d]*fbond_ji;
#endif
        }

        con.second[m_indexStretch] = s;
    }
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    // PD_bond
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);

    double dr_ij[m_dim];
    double energy = 0;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    for(auto &con:PDconnections)
    {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int col_j = m_pIds[id_j];

        const double vol_j = m_data(col_j, m_indexVolume);
        const double dr0         = con.second[m_indexDr0];
        const double volumeScaling   = con.second[m_indexVolumeScaling];
        const double g_ij = con.second[m_indexForceScaling];
        const double w_ij = con.second[m_indexWeightFunction];
        const double c_j = m_data(col_j, m_indexMicromodulus);
        const double c_ij = 0.5*(c_i + c_j)*w_ij*g_ij;

        double dr2 = 0;
        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(d, col_j) - m_r(d, i);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        energy += c_ij*(ds*ds)/dr0*vol_j*volumeScaling;
    }

    return 0.25*energy;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculatePotentialEnergy(const std::pair<int, int> &idCol, int indexPotential)
{
    // PD_bond
    int col_i = idCol.second;
    double vol_i = m_data(col_i, m_indexVolume);
    m_data(col_i, indexPotential) += calculatePotentialEnergyDensity(idCol)*vol_i;
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculateBondEnergy(const std::pair<int, int> &idCol, pair<int, vector<double> > &con)
{
    // PD_bond
    const int pId = idCol.first;
    const int col_i = idCol.second;
    const double c_i = m_data(col_i, m_indexMicromodulus);

    const int id_j = con.first;
    const int col_j = m_pIds[id_j];

    const double vol_j = m_data(col_j, m_indexVolume);
    const double dr0Len         = con.second[m_indexDr0];
    const double volumeScaling   = con.second[m_indexVolumeScaling];
    const double g_ij = con.second[m_indexForceScaling];
    const double w_ij = con.second[m_indexWeightFunction];
    const double c_j = m_data(col_j, m_indexMicromodulus);
    const double c = 0.5*(c_i + c_j)*w_ij*g_ij;

    double dr_ij[m_dim];
    double dr2 = 0;

    for(int d=0; d<m_dim; d++)
    {
        dr_ij[d] = m_r(d, col_j) - m_r(d, col_i);
        dr2 += dr_ij[d]*dr_ij[d];
    }

    const double drLen = sqrt(dr2);
    double ds = drLen - dr0Len;

    // To avoid roundoff errors
    if (fabs(ds) < THRESHOLD)
        ds = 0.0;

    const double energy = 0.5*c*(ds*ds)/dr0Len*vol_j*volumeScaling;

    return energy;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    // PD_bond
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    double r_i[m_dim];
    for(int d=0; d<m_dim; d++)
    {
        r_i[d] = m_r(d, i);
    }
    double dr_ij[m_dim];
    double f[m_dim];

    for(auto &con:PDconnections)
    {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double g_ij = con.second[m_indexForceScaling];
        const double w_ij = con.second[m_indexWeightFunction];
        const double c_ij = 0.5*(c_i + c_j)*w_ij*g_ij;
        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(d, j) - r_i[d];
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        double s = ds/dr0;
        //        double stretch = con.second[m_indexStretch];
        const double fbond_ij = c_ij*s*vol_j*volumeScaling/dr;
#ifdef USE_N3L
        const double bond_ji = c_ij*s*vol_i*volumeScaling/dr;
#endif
        for(int d=0; d<m_dim; d++)
        {
            f[d] = dr_ij[d]*fbond_ij;
        }
        m_data(i, indexStress[0]) += 0.5*f[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*f[Y]*dr_ij[Y];
        m_data(i, indexStress[3]) += 0.5*f[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[2]) += 0.5*f[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*f[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*f[Z]*dr_ij[Z];
        }
#ifdef USE_N3L
#ifdef USE_OPENMP
        //#pragma omp critical
#endif
        {
            m_data(j, indexStress[0]) += 0.5*bond_ji*dr_ij[X]*dr_ij[X];
            m_data(j, indexStress[1]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Y];
            m_data(j, indexStress[3]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Y];

            if(m_dim == 3)
            {
                m_data(j, indexStress[2]) += 0.5*bond_ji*dr_ij[Z]*dr_ij[Z];
                m_data(j, indexStress[4]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Z];
                m_data(j, indexStress[5]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Z];
            }
        }
#endif
    }
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    // PD_bond
    dt *= 1.1;
    const int pId = idCol.first;
    const int col_a = idCol.second;
    const double c_a = m_data(col_a, m_indexMicromodulus);

    double m[m_dim];
    double dR0[m_dim];

    const arma::mat & matR0 = m_particles.r0();
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
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_b = con.first;
            const int col_b = m_pIds[id_b];

            double sum = 0;
            for(int d=0; d<m_dim; d++)
            {
                dR0[d] = matR0(d, col_a) - matR0(d, col_b);
                sum += fabs(dR0[d]);
            }

            const double c_b = m_data(col_b, m_indexMicromodulus);
            const double vol_b = m_data(col_b, m_indexVolume);
            const double dr0 = con.second[m_indexDr0];
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double g_ab = con.second[m_indexForceScaling];
            const double w_ij = con.second[m_indexWeightFunction];
            const double c_ab = 0.5*(c_a + c_b)*w_ij*g_ab;
            const double dr0_3 = pow(dr0, 3);

            const double coeff = c_ab*volumeScaling*vol_b/dr0_3;
            k[i] += coeff*fabs(dR0[i])*sum;
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

    return 2*0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initialize(double E, double nu, double delta,
                                      int dim, double h)
{
    Force::initialize(E, nu, delta, dim, h);

    if(m_weightType == "constant")
    {
        initializeConstant();
    }
    else if(m_weightType == "linear")
    {
        initializeLinear();
    }
    else if(m_weightType == "gaussian")
    {
        initializeGaussian();
    }
    else if(m_weightType == "sigmoid")
    {
        initializeSigmoid();
    }
    else
    {
        std::cerr << "Weightfunction " << m_weightType << endl
                  << "for bondforce" << std::endl;
        exit(EXIT_FAILURE);
    }

}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeConstant()
{
    double c, k;

    if(m_dim == 3)
    {
        m_nu = 1./4.;
        k = m_E/(3.*(1. - 2.*m_nu));
        c = 18.*k/(M_PI*pow(m_delta, 4));
    }
    else if(m_dim == 2)
    {
        m_nu = 1./3.;
        k = m_E/(2.*(1. - m_nu));
        c = 12*k/(m_h*M_PI*pow(m_delta, 3));
    }
    else if(m_dim == 1)
    {
        m_nu = 1./4.;
        k = m_E;
        c = 2*m_E/(m_h*m_h*pow(m_delta, 2));
    }
    else
    {
        cerr << "ERROR: dimension " << m_dim << " not supported" << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    m_particles.setParameter("micromodulus", c);
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeLinear()
{
    double c, k;

    if(m_dim == 3)
    {
        m_nu = 1./4.;
        k = m_E/(3.*(1. - 2.*m_nu));
        c = 90.*k/(M_PI*pow(m_delta, 4));
    }
    else if(m_dim == 2)
    {
        m_nu = 1./3.;
        k = m_E/(2.*(1. - m_nu));
        c = 48.*k/(m_h*M_PI*pow(m_delta, 3));
    }
    else if(m_dim == 1)
    {
        cerr << "dim == 1 is not implemented for linear weightfunction" << endl;
        exit(EXIT_FAILURE);
    }

    m_particles.setParameter("micromodulus", c);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles.nParticles(); i++)
    {
        const int pId = i;
        const int col_i = i;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            const double dr0 = con.second[m_indexDr0];
            con.second[m_indexWeightFunction] = 1 - dr0/m_delta;
        }
    }
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeGaussian()
{
    double c, k;
    m_l = 0.5*m_delta;

    if(m_dim == 3)
    {
        m_nu = 1./4.;
        k = m_E/(3.*(1. - 2.*m_nu));
        c = 9.*k/(4*M_PI*pow(m_delta, 4)*(3. - 8.*exp(-1)));
    }
    else if(m_dim == 2)
    {
        m_nu = 1./3.;
        k = m_E/(2.*(1. - m_nu));
        c = 4.*k/(m_h*M_PI*pow(m_delta, 3)*(2. - 5.*exp(-1)));
    }
    else if(m_dim == 1)
    {
        m_nu = 1./4.;
        k = m_E;
        c = 2*k/(m_h*m_h*pow(m_delta, 2))*(1. - 2.*exp(-1));
        c = 2*k/(m_h*m_h*pow(m_delta, 2));
    }

    m_particles.setParameter("micromodulus", c);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles.nParticles(); i++)
    {
        const int pId = i;
        const int col_i = i;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            const double dr0Len = con.second[m_indexDr0];
            const double e_drl = exp(-dr0Len/m_l);

            con.second[m_indexWeightFunction] = e_drl;
        }
    }
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initializeSigmoid()
{
    double c, k;
    double alpha = 0.5*m_delta;
    double beta = m_delta/10.;

    if(m_dim == 3)
    {
        m_nu = 1./4.;
        k = m_E/(3.*(1. - 2.*m_nu));
        c = 18.*k/(M_PI*pow(m_delta, 4));
    }
    else if(m_dim == 2)
    {
        m_nu = 1./3.;
        k = m_E/(2.*(1. - m_nu));
        c = 12*k/(m_h*M_PI*pow(m_delta, 3));
    }
    else if(m_dim == 1)
    {
        m_nu = 1./4.;
        k = m_E;
        c = 2*m_E/(m_h*m_h*pow(m_delta, 2));
    }
    else
    {
        cerr << "ERROR: dimension " << m_dim << " not supported" << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    m_particles.setParameter("micromodulus", c);


#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles.nParticles(); i++)
    {
        const int pId = i;
        const int col_i = i;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            const double dr0 = con.second[m_indexDr0];
            const double e_drl = 1./(1. + exp((dr0 - alpha)/beta));

            con.second[m_indexWeightFunction] = e_drl;
        }
    }
}
//------------------------------------------------------------------------------
}
