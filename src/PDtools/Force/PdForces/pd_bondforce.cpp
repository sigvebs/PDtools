#include "pd_bondforce.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_bondForce::PD_bondForce(PD_Particles &particles):
    Force(particles, "PD bond force")
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexVolume = m_particles.getParamId("volume");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexForceScaling = m_particles.registerPdParameter("forceScalingBond", 1);
    m_indexStretch = m_particles.registerPdParameter("stretch");
    m_indexConnected = m_particles.getPdParamId("connected");
#if USE_N3L
    m_indexMyPdPosition = m_particles.getPdParamId("myPosistion");
#endif
    m_hasSurfaceCorrection = true;
    m_ghostParameters = {"volume", "micromodulus"};
    m_initialGhostParameters = {"volume", "micromodulus"};

    m_indexStress[0] = m_particles.registerParameter("s_xx");
    m_indexStress[1] = m_particles.registerParameter("s_yy");
    m_indexStress[2] = m_particles.registerParameter("s_xy");
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateForces(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
    const int nParticles = m_particles.nParticles();
#endif
    vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections_i.size();
    double dr_ij[m_dim];

    //----------------------------------
    // TMP - standard stres calc from
    m_data(i, m_indexStress[0]) = 0;
    m_data(i, m_indexStress[1]) = 0;
    m_data(i, m_indexStress[2]) = 0;
    //----------------------------------

    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con_i = PDconnections_i[l_j];
        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con_i.first;
        const int j = m_idToCol.at(id_j);

#if USE_N3L
        if(j<i)
            continue;
#endif
        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con_i.second[m_indexDr0];
        const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
        const double g_ij = con_i.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;
        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double ds = dr - dr0;
        const double s = ds/dr0;
        const double fbond_ij = c_ij*s*vol_j*volumeScaling_ij/dr;

        for(int d=0; d<m_dim; d++) {
            m_F(i, d) += dr_ij[d]*fbond_ij;
        }

        //----------------------------------
        // TMP - standard stres calc from
        m_data(i, m_indexStress[0]) += 0.5*dr_ij[0]*dr_ij[0]*fbond_ij;
        m_data(i, m_indexStress[1]) += 0.5*dr_ij[1]*dr_ij[1]*fbond_ij;
        m_data(i, m_indexStress[2]) += 0.5*dr_ij[0]*dr_ij[1]*fbond_ij;
        //----------------------------------

        con_i.second[m_indexStretch] = s;
#if USE_N3L
        if(j > i && j < nParticles) {
            const int myPos_j = con_i.second[m_indexMyPdPosition];
            vector<pair<int, vector<double>>> & PDconnections_j = m_particles.pdConnections(id_j);
            auto &con_j = PDconnections_j[myPos_j];
            const double volumeScaling_ji = con_j.second[m_indexVolumeScaling];
            const double fbond_ji = -c_ij*s*vol_i*volumeScaling_ji/dr;
            for(int d=0; d<m_dim; d++) {
                m_F(j, d) += dr_ij[d]*fbond_ji;
            }
            con_j.second[m_indexStretch] = s;
        }
#endif
    }
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateLinearForces(const int id_i, const int i)
{
    (void) id_i;
    (void) i;
}
//------------------------------------------------------------------------------
double PD_bondForce::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);
    double dr_ij[m_dim];
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    double energy = 0;
    for(auto &con:PDconnections) {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con.second[m_indexDr0];
        const double volumeScaling = con.second[m_indexVolumeScaling];
        const double c_j = m_data(j, m_indexMicromodulus);
        const double g_ij = con.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++) {
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
void PD_bondForce::calculatePotentialEnergy(const int id_i, const int i, int indexPotential)
{
    const double vol_i = m_data(i, m_indexVolume);
    m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i)*vol_i;
}
//------------------------------------------------------------------------------
void PD_bondForce::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{    
    const double c_i = m_data(i, m_indexMicromodulus);
#if USE_N3L
    const double vol_i = m_data(i, m_indexVolume);
    const int nParticles = m_particles.nParticles();
#endif

    const vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id_i);

    double dr_ij[m_dim];
    const int nConnections = PDconnections_i.size();

    for(int l_j=0; l_j<nConnections; l_j++) {
        const auto &con_i = PDconnections_i[l_j];
        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con_i.first;
        const int j = m_idToCol.at(id_j);
#if USE_N3L // Already computed
        if(j<i)
            continue;
#endif
        const double c_j = m_data(j, m_indexMicromodulus);
        const double vol_j = m_data(j, m_indexVolume);
        const double dr0 = con_i.second[m_indexDr0];
        const double volumeScaling = con_i.second[m_indexVolumeScaling];
        const double g_ij = con_i.second[m_indexForceScaling];
        const double c_ij = 0.5*(c_i + c_j)*g_ij;

        double dr2 = 0;

        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
        }

        const double dr = sqrt(dr2);
        const double s = (dr - dr0)/dr0;
        const double bond_ij = c_ij*s*vol_j*volumeScaling/dr;

        m_data(i, indexStress[0]) += 0.5*bond_ij*dr_ij[X]*dr_ij[X];
        m_data(i, indexStress[1]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Y];
        m_data(i, indexStress[2]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Y];

        if(m_dim == 3) {
            m_data(i, indexStress[3]) += 0.5*bond_ij*dr_ij[Z]*dr_ij[Z];
            m_data(i, indexStress[4]) += 0.5*bond_ij*dr_ij[X]*dr_ij[Z];
            m_data(i, indexStress[5]) += 0.5*bond_ij*dr_ij[Y]*dr_ij[Z];
        }
#if USE_N3L
        if(j > i && j < nParticles) {
            const int myPos_j = con_i.second[m_indexMyPdPosition];
            const vector<pair<int, vector<double>>> & PDconnections_j = m_particles.pdConnections(id_j);
            const auto &con_j = PDconnections_j[myPos_j];
            const double volumeScaling_ji = con_j.second[m_indexVolumeScaling];
            const double bond_ji = c_ij*s*vol_i*volumeScaling_ji/dr;
            m_data(j, indexStress[0]) += 0.5*bond_ji*dr_ij[X]*dr_ij[X];
            m_data(j, indexStress[1]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Y];
            m_data(j, indexStress[2]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Y];

            if(m_dim == 3) {
                m_data(j, indexStress[3]) += 0.5*bond_ji*dr_ij[Z]*dr_ij[Z];
                m_data(j, indexStress[4]) += 0.5*bond_ji*dr_ij[X]*dr_ij[Z];
                m_data(j, indexStress[5]) += 0.5*bond_ji*dr_ij[Y]*dr_ij[Z];
            }
        }
#endif
    }
}
//------------------------------------------------------------------------------
double PD_bondForce::calculateStableMass(const int id_a, const int a, double dt)
{
    dt *= 1.1;

    const arma::mat & R0 = m_particles.r0();
    const double c_a = m_data(a, m_indexMicromodulus);

    double m[m_dim];
    double dr0[m_dim];

    for(int d=0; d<m_dim;d++) {
        m[d] = 0;
    }

    const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_a);

    double k[m_dim];

    for(int i=0; i<m_dim; i++) {
        for(int d=0; d<m_dim;d++) {
            k[d] = 0;
        }
        for(auto &con:PDconnections)
        {
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_b = con.first;
            const int b = m_idToCol.at(id_b);

            for(int d=0; d<m_dim;d++) {
                dr0[d] =   R0(a, d) - R0(b, d);
            }

            const double c_b = m_data(b, m_indexMicromodulus);
            const double vol_b = m_data(b, m_indexVolume);
            const double dr0Len = con.second[m_indexDr0];
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double g_ab = con.second[m_indexForceScaling];
            const double c_ab = 0.5*(c_a + c_b)*g_ab;

            const double dr0Len3 = pow(dr0Len, 3);
            const double coeff = c_ab*volumeScaling*vol_b/dr0Len3;

            double sum = 0;

            for(int d=0; d<m_dim; d++) {
                sum += fabs(dr0[d]);
            }
            sum *= fabs(dr0[i])*coeff;

            k[i] += sum;
        }

        m[i] = k[i];
    }

    double stiffness = 0;

    for(int d=0;d<m_dim; d++) {
        if(m[d]>stiffness) {
            stiffness = m[d];
        }
    }

    return 4.*0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void PD_bondForce::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    double k;
    double c;

    if(dim == 3) {
        nu = 1./4.;
        k = E/(3.*(1. - 2.*nu));
        c = 18.*k/(M_PI*pow(delta, 4));
    } else if(dim == 2) {
        nu = 1./3.;
        k = E/(2.*(1. - nu));
        c = 12*k/(h*M_PI*pow(delta, 3));
    } else if(dim == 1) {
        nu = 1./4.;
        k = E;
        c = 2*E/(h*h*pow(delta, 2));
    } else {
        cerr << "ERROR: dimension " << dim << " not supported" << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    if(m_numericalInitialization){
        m_particles.setParameter("micromodulus", k);
    } else {
        m_particles.setParameter("micromodulus", c);
        return;
    }

    /*
    // Scaling the bond force
    const double dimScaling = 2.*k*pow(dim, 2.);
    const ivec & colToId = m_particles.colToId();
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++) {
        const int pId = colToId(i);
        double dRvolume = 0;
        double v = 0;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int col_j = m_idToCol.at(id_j);
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double dr0Len = con.second[m_indexDr0];
            const double g_ij = con.second[m_indexForceScaling];
            dRvolume +=   dr0Len*m_data(col_j, m_indexVolume)*volumeScaling*g_ij;
            v += m_data(col_j, m_indexVolume)*volumeScaling;
        }
        const double c_i = dimScaling/dRvolume;
        m_data(i, m_indexMicromodulus) = c_i;
//        double l_delta = pow((3.*v/(4.*M_PI)), 1./3.);
//        double c_j = 18*k/(M_PI*pow(l_delta, 4));
//        m_data(col_i, m_indexMicromodulus) = c_j;
    }
    */
}
//------------------------------------------------------------------------------
}
