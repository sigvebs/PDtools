#include "pd_bondforcegaussian.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondforceGaussian::PD_bondforceGaussian(PD_Particles &particles):
    Force(particles),
    m_r(m_particles.r()),
    m_F(m_particles.F()),
    m_data(m_particles.data()),
    m_pIds(m_particles.pIds())
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexExponent = m_particles.registerPdParameter("exponent");
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
    m_indexStretch = m_particles.registerPdParameter("stretch");
}
//------------------------------------------------------------------------------
PD_bondforceGaussian::~PD_bondforceGaussian()
{

}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::calculateForces(const std::pair<int, int> &idCol)
{
    // PD_bond
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double x_i = m_r(X, col_i);
    double y_i = m_r(Y, col_i);
    double z_i = m_r(Z, col_i);
    double dr_ij[3];

    for(auto &con:PDconnections)
    {
        int id_j = con.first;
        int col_j = m_pIds[id_j];

        double c_j = m_data(col_j, m_indexMicromodulus);
        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len         = con.second[m_indexDr0];
        double volumeScaling   = con.second[m_indexVolumeScaling];
        double exponential = con.second[m_indexExponent];
        double c_ij = 0.5*(c_i + c_j)*exponential;

        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        double stretch = ds/dr0Len;
        double fbond = c_ij*stretch*vol_j*volumeScaling/drLen;

        m_F(X, col_i) += dr_ij[X]*fbond;
        m_F(Y, col_i) += dr_ij[Y]*fbond;
        m_F(Z, col_i) += dr_ij[Z]*fbond;

        con.second[m_indexStretch] = stretch;
    }
}
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    // PD_bond
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    double dr_ij[3];
    double x_i = m_r(X, col_i);
    double y_i = m_r(Y, col_i);
    double z_i = m_r(Z, col_i);

    double energy = 0;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    for(auto &con:PDconnections)
    {

        int id_j = con.first;
        int col_j = m_pIds[id_j];

        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len         = con.second[m_indexDr0];
        double volumeScaling   = con.second[m_indexVolumeScaling];
        double exponential = con.second[m_indexExponent];
        double c_j = m_data(col_j, m_indexMicromodulus);

        double c = 0.5*(c_i + c_j)*exponential;

        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        energy += c*(ds*ds)/dr0Len*vol_j*volumeScaling;
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
void PD_bondforceGaussian::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    // PD_bond
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double x_i = m_r(X, col_i);
    double y_i = m_r(Y, col_i);
    double z_i = m_r(Z, col_i);
    double dr_ij[3];
    double f[3];

    for(auto &con:PDconnections)
    {
        int id_j = con.first;
        int col_j = m_pIds[id_j];

        double c_j = m_data(col_j, m_indexMicromodulus);
        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len = con.second[m_indexDr0];
        double volumeScaling = con.second[m_indexVolumeScaling];
        double exponential = con.second[m_indexExponent];
        double c_ij = 0.5*(c_i + c_j)*exponential;

        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        double stretch = ds/dr0Len;
        stretch = con.second[m_indexStretch];
        double fbond = c_ij*stretch*vol_j*volumeScaling/drLen;

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
//------------------------------------------------------------------------------
double PD_bondforceGaussian::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    // PD_bond
    dt *= 1.1;
    int pId = idCol.first;
    int col_a = idCol.second;
    double c_a = m_data(col_a, m_indexMicromodulus);

    double m[3];
    double dR0[3];
    m[X] = 0;
    m[Y] = 0;
    m[Z] = 0;

    const arma::mat & matR0 = m_particles.r0();

    double x_a = matR0(X, col_a);
    double y_a = matR0(Y, col_a);
    double z_a = matR0(Z, col_a);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double k[3];

    for(int i=0; i<3; i++)
    {
        k[X] = 0;
        k[Y] = 0;
        k[Z] = 0;


        for(auto &con:PDconnections)
        {
            int id_b = con.first;
            int col_b = m_pIds[id_b];

            dR0[X] = x_a - matR0(X, col_b);
            dR0[Y] = y_a - matR0(Y, col_b);
            dR0[Z] = z_a - matR0(Z, col_b);

            double c_b = m_data(col_b, m_indexMicromodulus);
            double vol_b = m_data(col_b, m_indexVolume);
            double dr0Len = con.second[m_indexDr0];
            double volumeScaling = con.second[m_indexVolumeScaling];
            double exponential = con.second[m_indexExponent];
            double c_ab = 0.5*(c_a + c_b)*exponential;
            double dr0Len2 = pow(dr0Len, 2);

            k[i] += c_ab*volumeScaling*fabs(dR0[i])*vol_b/dr0Len2;
        }

        m[i] = k[i];
    }

    double stiffness = 0;

    for(int d=0;d<3; d++)
    {
        if(m[d]>stiffness)
        {
            stiffness = m[d];
        }
    }

    return 0.5*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void PD_bondforceGaussian::initialize(double E, double nu, double delta, int dim, double h)
{
    double k;
    m_l = 0.5*delta;

    if(dim == 3)
    {
        nu = 1./4.;
        k = E/(3.*(1. - 2.*nu));
    }
    else if(dim == 2)
    {
        nu = 1./3.;
        k = E/(2.*(1. - nu));
    }

    double dimScaling = 2.*pow(dim, 2.)*k;

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(int i=0; i<m_particles.nParticles(); i++)
    {
        int pId = i;
        int col_i = i;
        double dRvolume = 0;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];
            double volumeScaling = con.second[m_indexVolumeScaling];
            double volume = con.second[m_indexVolume];
            double dr0Len = con.second[m_indexDr0];
            double e_drl = exp(-dr0Len/m_l);

            dRvolume += dr0Len*volume;
            //            dRvolume += dr0Len*volume*volumeScaling;
            con.second[m_indexExponent] = e_drl;
        }
        double c = dimScaling/dRvolume;
        m_data(col_i, m_indexMicromodulus) = c;
    }
}
//------------------------------------------------------------------------------
}
