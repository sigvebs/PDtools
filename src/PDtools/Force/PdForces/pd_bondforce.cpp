#include "pd_bondforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondForce::PD_bondForce(PD_Particles &particles):
    Force(particles)
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
    m_indexForceScaling = m_particles.registerPdParameter("forceScalingBond", 1);
    //m_indexWeightfunction = m_particles.registerPdParameter("weightFunction", 1);
    m_indexStretch = m_particles.registerPdParameter("stretch");
    m_indexConnected = m_particles.getPdParamId("connected");
    m_indexCompute = m_particles.getPdParamId("compute");
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

    const arma::vec3 r_i = m_r.row(col_i);

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

        const double dx = m_r(col_j, X) - r_i(X);
        const double dy = m_r(col_j, Y) - r_i(Y);
        const double dz = m_r(col_j, Z) - r_i(Z);

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
    m_F(col_i, X) += f_ix;
    m_F(col_i, Y) += f_iy;
    m_F(col_i, Z) += f_iz;
#else
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);
#ifdef USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
#endif

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5 || !con.second[m_indexCompute])
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0         = con.second[m_indexDr0];
        const double volumeScaling   = con.second[m_indexVolumeScaling];
        const double g_ij = con.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        double dr = sqrt(dr2);
        double ds = dr - dr0;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        if(dr < 0)
        {
            ds = 0;
            dr = 1.0;
        }
        const double s = ds/dr0;
        const double fbond_ij = c_ij*s*vol_j*volumeScaling/dr;
#ifdef USE_N3L
        const double fbond_ji = -c_ij*s*vol_i*volumeScaling/dr;
#endif

        for(int d=0; d<m_dim; d++)
        {
            m_F(i, d) += dr_ij[d]*fbond_ij;
#ifdef USE_N3L
            m_F(j, d) += dr_ij[d]*fbond_ji;
#endif
        }

        con.second[m_indexStretch] = s;
    }
#endif 
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateLinearForces(const std::pair<int, int> &idCol)
{
    (void) idCol;
}
//------------------------------------------------------------------------------
double PD_bondForce::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);
    double dr_ij[m_dim];

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double energy = 0;
    for(auto &con:PDconnections)
    {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double vol_j = m_data(j, m_indexVolume);
        const double dr0         = con.second[m_indexDr0];
        const double volumeScaling   = con.second[m_indexVolumeScaling];
        const double c_j = m_data(j, m_indexMicromodulus);
        const double g_ij = con.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double s = (dr - dr0)/dr0;
        energy += c_ij*(s*s)*dr0*vol_j*volumeScaling;
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
    const int pId = idCol.first;
    const int i = idCol.second;
    const double c_i = m_data(i, m_indexMicromodulus);
#ifdef USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
#endif

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double dr_ij[m_dim];
    const int nConnections = PDconnections.size();

    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_indexConnected] <= 0.5 || !con.second[m_indexCompute])
            continue;

        const int id_j = con.first;
        const int j = m_pIds[id_j];

        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        //const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double g_ij = con.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++)
        {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
//        const double s = (dr - dr0)/dr0;
        const double s = con.second[m_indexStretch];
        const double bond_ij = c_ij*s*vol_j*volumeScaling/dr;
#ifdef USE_N3L
        const double bond_ji = c_ij*s*vol_i*volumeScaling/dr;
#endif
        m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
        m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

        if(m_dim == 3)
        {
            m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
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
double PD_bondForce::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    dt *= 1.1;
    const int pId = idCol.first;
    const int col_a = idCol.second;
    const double c_a = m_data(col_a, m_indexMicromodulus);

    double m[3];
    double dr0[3];
    m[X] = 0;
    m[Y] = 0;
    m[Z] = 0;

    const arma::mat & matR0 = m_particles.r0();

    const double x_a = matR0(col_a, X);
    const double y_a = matR0(col_a, Y);
    const double z_a = matR0(col_a, Z);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

    double k[3];

    for(int i=0; i<3; i++)
    {
        k[X] = 0;
        k[Y] = 0;
        k[Z] = 0;

        for(auto &con:PDconnections)
        {
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_b = con.first;
            const int col_b = m_pIds[id_b];

            dr0[X] = x_a - matR0(col_b, X);
            dr0[Y] = y_a - matR0(col_b, Y);
            dr0[Z] = z_a - matR0(col_b, Z);

            const double c_b = m_data(col_b, m_indexMicromodulus);
            const double vol_b = m_data(col_b, m_indexVolume);
            const double dr0Len = con.second[m_indexDr0];
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double g_ab = con.second[m_indexForceScaling];
            const double c_ab = 0.5*(c_a + c_b)*g_ab;

            const double dr0Len3 = pow(dr0Len, 3);
            const double coeff = c_ab*volumeScaling*vol_b/dr0Len3;

            double sum = 0;

            for(int j=0; j<3; j++)
            {
                sum += fabs(dr0[j]);
            }
            sum *= fabs(dr0[i])*coeff;

            k[i] += sum;
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

    double stableMass = 4*0.25*pow(dt, 2)*stiffness;
    return stableMass;
}
//------------------------------------------------------------------------------
void PD_bondForce::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
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
    else if(dim == 1)
    {
        nu = 1./4.;
        k = E;
        c = 2*E/(h*h*pow(delta, 2));
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported" << endl;
        cerr << "use 1, 2 or 3." << endl;
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

    // Scaling the bond force
    const double dimScaling = 2.*k*pow(dim, 2.);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        const int pId = i;
        const int col_i = i;
        double dRvolume = 0;
        double v = 0;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = m_pIds[id_j];
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double dr0Len = con.second[m_indexDr0];
            const double g_ij = con.second[m_indexForceScaling];
            dRvolume +=   dr0Len*m_data(col_j, m_indexVolume)*volumeScaling*g_ij;
            v += m_data(col_j, m_indexVolume)*volumeScaling;
        }
        const double c_i = dimScaling/dRvolume;
        m_data(col_i, m_indexMicromodulus) = c_i;
//        double l_delta = pow((3.*v/(4.*M_PI)), 1./3.);
//        double c_j = 18*k/(M_PI*pow(l_delta, 4));
//        m_data(col_i, m_indexMicromodulus) = c_j;
    }
}
//------------------------------------------------------------------------------
}
