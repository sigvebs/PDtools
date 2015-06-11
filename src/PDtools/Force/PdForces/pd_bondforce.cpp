#include "pd_bondforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondForce::PD_bondForce(PD_Particles &particles):
    Force(particles),
    m_r(m_particles.r()),
    m_F(m_particles.F()),
    m_data(m_particles.data()),
    m_pIds(m_particles.pIds())
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);

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
PD_bondForce::~PD_bondForce()
{
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateForces(const pair<int, int> &idCol)
{
#if 0
    const int pId = idCol.first;
    const int col_i = idCol.second;
    const double c_i = m_data(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    const arma::vec3 r_i = m_r.col(col_i);

    double f_ix = 0;
    double f_iy = 0;
    double f_iz = 0;

    const int nConnections = PDconnections.size();
    double col_js[nConnections];
    double vol_js[nConnections];
    double c_ijs[nConnections];
    double dr0Lens[nConnections];
    double volumeScalings[nConnections];
    double stretchs[nConnections];

    for(int j=0; j<nConnections; j++)
    {
        const auto &con = PDconnections.at(j);
        const int id_j = con.first;
        const int col_j = m_pIds.at(id_j);
        col_js[j] = col_j;
        vol_js[j] = m_data(col_j, m_indexVolume);
        c_ijs[j] = 0.5*(c_i + m_data(col_j, m_indexMicromodulus));
        dr0Lens[j] = con.second[m_indexDr0];
        volumeScalings[j] = con.second[m_indexVolumeScaling];
    }

#pragma simd reduction(+: f_ix, f_iy, f_iz)
    for(int j=0; j<nConnections; j++)
    {
        const int col_j = col_js[j];
        const double dr0Len = dr0Lens[j];

        const double dx = m_r(X, col_j) - r_i(X);
        const double dy = m_r(Y, col_j) - r_i(Y);
        const double dz = m_r(Z, col_j) - r_i(Z);

        const double drLen = sqrt(dx*dx + dy*dy + dz*dz);
        const double ds = drLen - dr0Len;

        double stretch = ds/dr0Len;
        stretch = fabs(ds) < THRESHOLD ? 0:ds;
        const double fbond = c_ijs[j]*stretch*vol_js[j]*volumeScalings[j]/drLen;

        f_ix += dx*fbond;
        f_iy += dy*fbond;
        f_iz += dz*fbond;
        stretchs[j] = stretch;
    }
    for(int j=0; j<nConnections; j++)
    {
        PDconnections[j].second[m_indexStretch] = stretchs[j];
    }
    m_F(X, col_i) += f_ix;
    m_F(Y, col_i) += f_iy;
    m_F(Z, col_i) += f_iz;
#else
    int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double r_i[m_dim];
    double f_i[m_dim];
    for(int d=0; d<m_dim; d++)
    {
        f_i[d] = 0;
        r_i[d] = m_r(d, col_i);
    }

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    for(int j=0; j<nConnections; j++)
    {
        auto &con = PDconnections[j];
        int id_j = con.first;
//        int col_j = col_js[j];
        int col_j = m_pIds[id_j];

        double c_j = m_data(col_j, m_indexMicromodulus);
        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len         = con.second[m_indexDr0];
        double volumeScaling   = con.second[m_indexVolumeScaling];
        double c_ij = 0.5*(c_i + c_j);

        double drSquared = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(d, col_j) - r_i[d];
            drSquared += dr_ij[d]*dr_ij[d];
        }

        double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        double stretch = ds/dr0Len;
        double fbond = c_ij*stretch*vol_j*volumeScaling/drLen;

        for(int d=0; d<m_dim; d++)
        {
            f_i[d] += dr_ij[d]*fbond;
        }

        con.second[m_indexStretch] = stretch;
    }

    for(int d=0; d<m_dim; d++)
    {
        m_F(d, col_i) += f_i[d];
    }
#endif
    /*
     *     int pId = idCol.first;
    int col_i = idCol.second;
    double c_i = m_data(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

//    double x_i = m_r(X, col_i);
//    double y_i = m_r(Y, col_i);
//    double z_i = m_r(Z, col_i);
    double r_i[3];
    for(int d=0; d<3; d++)
    {
        r_i[d] = m_r(d, col_i);
    }
    double dr_ij[3];

    const int nConnections = PDconnections.size();
#pragma simd
    for(int j=0; j<nConnections; j++)
//        for(auto &con:PDconnections)
    {
        auto &con = PDconnections[j];
        int id_j = con.first;
        int col_j = m_pIds[id_j];

        double c_j = m_data(col_j, m_indexMicromodulus);
        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len         = con.second[m_indexDr0];
        double volumeScaling   = con.second[m_indexVolumeScaling];
        double c_ij = 0.5*(c_i + c_j);

        double drSquared = 0;
        for(int d=0; d<3; d++)
        {
            dr_ij[d] = m_r(d, col_j) - r_i[d];
            drSquared = dr_ij[d]*dr_ij[d];
        }
//        dr_ij[X] = m_r(X, col_j) - x_i;
//        dr_ij[Y] = m_r(Y, col_j) - y_i;
//        dr_ij[Z] = m_r(Z, col_j) - z_i;

//        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
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
    }*/
}
//------------------------------------------------------------------------------
double PD_bondForce::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
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
        double c_j = m_data(col_j, m_indexMicromodulus);

        double c = 0.5*(c_i + c_j);

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
void PD_bondForce::calculatePotentialEnergy(const std::pair<int, int> &idCol, int indexPotential)
{
    int col_i = idCol.second;
    double vol_i = m_data(col_i, m_indexVolume);
    m_data(col_i, indexPotential) += calculatePotentialEnergyDensity(idCol)*vol_i;
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
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
        double c_ij = 0.5*(c_i + c_j);

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
double PD_bondForce::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    dt *= 1.1;
#if 0
    if(0)
    {

        int pId = idCol.first;
        int col_i = idCol.second;
        double c_i = m_data(col_i, m_indexMicromodulus);
        double stiffness_xyz[3];
        double dR0[3];
        stiffness_xyz[X] = 0;
        stiffness_xyz[Y] = 0;
        stiffness_xyz[Z] = 0;

        const arma::mat & matR0 = m_particles.r0();

        double x_i = matR0(X, col_i);
        double y_i = matR0(Y, col_i);
        double z_i = matR0(Z, col_i);

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);


        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];

            dR0[X] = x_i - matR0(X, col_j);
            dR0[Y] = y_i - matR0(Y, col_j);
            dR0[Z] = z_i - matR0(Z, col_j);

            double c_j = m_data(col_j, m_indexMicromodulus);
            double vol_j = m_data(col_j, m_indexVolume);
            double dr0Len = con.second[m_indexDr0];
            double volumeScaling = con.second[m_indexVolumeScaling];
            double c_ij = 0.5*(c_i + c_j);

            double dR0xyz = abs(dR0[X]) + abs(dR0[Y]) + abs(dR0[Z]);
            double dr0Len3 = pow(dr0Len, 3);
            for(int d=0; d<3; d++)
            {
                stiffness_xyz[d] +=
                        c_ij
                        * abs(dR0[d])
                        * dR0xyz
                        * vol_j
                        * volumeScaling
                        / dr0Len3;
            }
        }

        double stiffness = 0;

        for(int d=0;d<3; d++)
        {
            if(stiffness_xyz[d]>stiffness)
            {
                stiffness = stiffness_xyz[d];
            }
        }

        return 0.25*pow(dt, 2)*stiffness;
    }
    else
    {
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

        double k_i[3];

        for(int i=0; i<3; i++)
        {
            k_i[X] = 0;
            k_i[Y] = 0;
            k_i[Z] = 0;

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
                double c_ab = 0.5*(c_a + c_b);

                double dr0Len3 = pow(dr0Len, 3);

                for(int j=0; j<3; j++)
                {
                    k_i[j] +=
                            c_ab
                            * dR0[i]*dR0[j]
                            * vol_b
                            * volumeScaling
                            / dr0Len3;
                }
            }

            for(int j=0; j<3; j++)
            {
                m[i] += fabs(k_i[j]);
            }
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
//        return 0.25*pow(dt, 2)*stiffness;
    }
#else if
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

    double k_one[3];
    double k_two[3];


    for(int i=0; i<3; i++)
    {
        k_one[X] = 0;
        k_one[Y] = 0;
        k_one[Z] = 0;

        k_two[X] = 0;
        k_two[Y] = 0;
        k_two[Z] = 0;

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
            double c_ab = 0.5*(c_a + c_b);

            double dr0Len3 = pow(dr0Len, 3);
            double coeff = c_ab*volumeScaling*vol_b/dr0Len3;
            double term1 = 0;
            double term2 = 0;

            for(int j=0; j<3; j++)
            {
                term1 += dR0[i]*dR0[j];
                term2 += fabs(dR0[j]);
            }
            term1 *= coeff;
            term2 *= fabs(dR0[i])*coeff;

            k_one[i] += term1;
            k_two[i] += term2;
        }

//        m[i] = fabs(k_one[i]) - k_two[i];
        m[i] = k_two[i];
    }

    double stiffness = 0;

    for(int d=0;d<3; d++)
    {
        if(m[d]>stiffness)
        {
            stiffness = m[d];
        }
    }

    return 0.25*pow(dt, 2)*stiffness;
#endif
}
//------------------------------------------------------------------------------
void PD_bondForce::initialize(double E, double nu, double delta, int dim, double h)
{
    double k;
    double c;

    if(dim == 3)
    {
        nu = 1./4.;
        k = E/(3.*(1. - 2.*nu));
        c = 18.*k/(M_PI*pow(delta, 4));
    }
    else if(dim == 2)
    {
        nu = 1./3.;
        k = E/(2.*(1. - nu));
        c = 12*k/(h*M_PI*pow(delta, 3));
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported" << endl;
        cerr << "use 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    if(m_numericalInitialization){
        m_particles.setParameter("micromodulus", k);
    }
    else
    {
        m_particles.setParameter("micromodulus", c);
        return;
    }

    double dimScaling = 2.*pow(dim, 2.);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
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
            double dr0Len = con.second[m_indexDr0];

            dRvolume += dr0Len*m_data(col_j, m_indexVolume);
            //            dRvolume += dr0Len*data(col_j, indexVolume)*volumeScaling;
        }
        m_data(col_i, m_indexMicromodulus) *= dimScaling/dRvolume;
    }
}
//------------------------------------------------------------------------------
}

