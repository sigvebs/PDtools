#include "pd_pmb.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_PMB::PD_PMB(PD_Particles &particles, double lc, double delta, double alpha):
    Force(particles)
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexS0 = m_particles.getParamId("s0");
    m_indexS00 = m_particles.registerPdParameter("s00");
    m_indexConnected = m_particles.getPdParamId("connected");
    m_indexS_new = m_particles.registerParameter("s_new");
    m_indexUnbreakable = m_particles.registerParameter("unbreakable");


    m_indexMicromodulus = m_particles.registerParameter("micromodulus");

    m_lc = lc;
    m_delta = delta;
    m_alpha = alpha;
    m_indexVolume = m_particles.getParamId("volume");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexForceScaling = m_particles.registerPdParameter("forceScalingBond", 1.);
    m_indexStretch = m_particles.registerPdParameter("stretch");

    int nParticles = particles.nParticles();
    f = new double* [nParticles];
    x = new double* [nParticles];
    r0 = new double* [nParticles];

    for(int i=0;i<m_dim; i++) {
        f[i] = m_F.colptr(i);
        x[i] = m_r.colptr(i);
        r0[i] = m_r0.colptr(i);
    }

    m_ghostParameters.push_back("volume");
    m_ghostParameters.push_back("unbreakable");
    m_ghostParameters.push_back("micromodulus"); // For contact forces
    m_initialGhostParameters = {"volume", "s0"};
}
//------------------------------------------------------------------------------
PD_PMB::~PD_PMB()
{
}
//------------------------------------------------------------------------------
void PD_PMB::calculateForces(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
//    vector<pair<int, vector<double>> *> removeParticles;

    int jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq, r, dr, rk, fbond;
    double delta, stretch, s0_new;

    jnum = PDconnections.size();
    const double lc = m_lc;
    delta = m_delta;
    const double half_lc = 0.5*lc;
    double vfrac_scale = 1.0;

//    bool first;

//    xtmp = x[i][0];
//    ytmp = x[i][1];
//    ztmp = x[i][2];

    xtmp = m_r(i,0);
    ytmp = m_r(i,1);
    ztmp = m_r(i,2);

    m_data(i, m_indexS0) = m_data(i, m_indexS_new);
    const double s0_i = m_data(i, m_indexS0);
    s0_new = numeric_limits<double>::max();
//    first = true;

    for (int jj = 0; jj < jnum; jj++) {
        auto &con = PDconnections[jj];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);
        const double dr0 = con.second[m_indexDr0];

        // compute force density, add to PD equation of motion
        delx = xtmp - m_r(j,0);
        dely = ytmp - m_r(j,1);
        delz = ztmp - m_r(j,2);
//        delx = xtmp - x[j][0];
//        dely = ytmp - x[j][1];
//        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);
        dr = r - dr0;

        // avoid roundoff errors
        if (fabs(dr) < 2.2204e-016) {
            dr = 0.0;
        }

        // scale vfrac[j] if particle j near the horizon
        const double c_j = m_data(j, m_indexMicromodulus);
        const double c_ij = 0.5*(c_i + c_j);

        if ((fabs(dr0 - delta)) <= half_lc) {
            vfrac_scale = (-1.0/(2*half_lc))*(dr0)
                    + (1.0 + ((delta - half_lc)/(2*half_lc) ) );
        }
        else {
            vfrac_scale = 1.0;
        }

        const double vfrac_j = m_data(j, m_indexVolume);

        stretch = dr/dr0;
        rk = (c_ij * vfrac_j) * vfrac_scale * stretch;

        if (r > 0.0)
            fbond = -(rk/r);
        else
            fbond = 0.0;

        m_F(i,X) += delx*fbond;
        m_F(i,Y) += dely*fbond;
        m_F(i,Z) += delz*fbond;
//        f[i][0] += delx*fbond;
//        f[i][1] += dely*fbond;
//        f[i][2] += delz*fbond;

        const double s0_j = m_data(j, m_indexS0);

        if (stretch > min(s0_i, s0_j)){
            if( m_data(i, m_indexUnbreakable) > 0|| m_data(j, m_indexUnbreakable) > 0) {

            } else {
                cout << "broken:" << id_i << " " << id_j << " s:" << stretch
                     << " s0_i:" << s0_i  << " s0_j:" << s0_i << endl;
                con.second[m_indexConnected] = 0;
            }
        }
//            removeParticles.push_back(&con);

        // update s0 for next timestep
        const double s00 = con.second[m_indexS00];
        s0_new = s00;
        /*
        if (first)
           s0_new = s00 - (m_alpha * stretch);
        else
           s0_new = max(s0_new, s00 - m_alpha * stretch);
        first = false;
        */
      }

    m_data(i, m_indexS_new) = s0_new;
}
//------------------------------------------------------------------------------
double PD_PMB::calculatePotentialEnergyDensity(const int id_i, const int i)
{
    const double c_i = m_data(i, m_indexMicromodulus);

    double dr_ij[3];
    const double x_i = m_r(i, X);
    const double y_i = m_r(i, Y);
    const double z_i = m_r(i, Z);

    double energy = 0;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    for(auto &con:PDconnections) {
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_indexVolume);
        const double dr0Len         = con.second[m_indexDr0];
        const double volumeScaling   = con.second[m_indexVolumeScaling];
        const double c_j = m_data(j, m_indexMicromodulus);
        const double g_ij = con.second[m_indexForceScaling];

        const double c_ij = 0.5*(c_i + c_j)*g_ij;


        dr_ij[X] = m_r(j, X) - x_i;
        dr_ij[Y] = m_r(j, Y) - y_i;
        dr_ij[Z] = m_r(j, Z) - z_i;

        const double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        const double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        energy += c_ij*(ds*ds)/dr0Len*vol_j*volumeScaling;
    }

    return 0.25*energy;
}
//------------------------------------------------------------------------------
void PD_PMB::calculatePotentialEnergy(const int id_i, const int i, int indexPotential)
{
    const double vol_i = m_data(i, m_indexVolume);
    m_data(i, indexPotential) += calculatePotentialEnergyDensity(id_i, i)*vol_i;
}
//------------------------------------------------------------------------------
void PD_PMB::calculateStress(const int id_i, const int i, const int (&indexStress)[6])
{
    const double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    vector<pair<int, vector<double>> *> removeParticles;

    // loop over my particles and their partners
    // partner list contains all bond partners, so I-J appears twice
    // if bond already broken, skip this partner
    // first = true if this is first neighbor of particle i
    int jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq, r, dr, rk, fbond;
    double delta, stretch;
//    double s0_new;

    jnum = PDconnections.size();
    double lc = m_lc;
    delta = m_delta;
    const double half_lc = 0.5*lc;
    double vfrac_scale = 1.0;

//    bool first; 
    //    xtmp = x[i][0];
    //    ytmp = x[i][1];
    //    ztmp = x[i][2];

    xtmp = m_r(i,0);
    ytmp = m_r(i,1);
    ztmp = m_r(i,2);

    double f_tmp[3];

    m_data(i, m_indexS0) = m_data(i, m_indexS_new);
//    s0_new = numeric_limits<double>::max();
//    first = true;

    for (int jj = 0; jj < jnum; jj++) {
        auto &con = PDconnections[jj];
        if(con.second[m_indexConnected] <= 0.5)
            continue;
        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);
        const double dr0 = con.second[m_indexDr0];

        // compute force density, add to PD equation of motion
        delx = xtmp -  m_r(j,0);
        dely = ytmp -  m_r(j,1);
        delz = ztmp -  m_r(j,2);

        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);
        dr = r - dr0;

        // avoid roundoff errors
        if (fabs(dr) < 2.2204e-016) {
            dr = 0.0;
        }

        // scale vfrac[j] if particle j near the horizon
        const double c_j = m_data(j, m_indexMicromodulus);
        const double c_ij = 0.5*(c_i + c_j);

        if ((fabs(dr0 - delta)) <= half_lc) {
            vfrac_scale = (-1.0/(2*half_lc))*(dr0)
                    + (1.0 + ((delta - half_lc)/(2*half_lc) ) );
        } else {
            vfrac_scale = 1.0;
        }

        const double vfrac_j = m_data(j, m_indexVolume);

        stretch = dr / dr0;
        rk = (c_ij * vfrac_j) * vfrac_scale * stretch;

        if (r > 0.0)
            fbond = -(rk/r);
        else
            fbond = 0.0;

        f_tmp[X] = delx*fbond;
        f_tmp[Y] = dely*fbond;
        f_tmp[Z] = delz*fbond;

        m_data(i, indexStress[0]) += 0.5*f_tmp[X]*delx;
        m_data(i, indexStress[1]) += 0.5*f_tmp[Y]*dely;
        m_data(i, indexStress[2]) += 0.5*f_tmp[X]*dely;

        if(m_dim == 3) {
            m_data(i, indexStress[3]) += 0.5*f_tmp[Z]*delz;
            m_data(i, indexStress[4]) += 0.5*f_tmp[X]*delz;
            m_data(i, indexStress[5]) += 0.5*f_tmp[Y]*delz;
        }
      }
}
//------------------------------------------------------------------------------
double PD_PMB::calculateStableMass(const int id_a, const int a, double dt)
{
    dt *= 1.1;
    const double c_a = m_data(a, m_indexMicromodulus);

    double m[3];
    double dR0[3];
    m[X] = 0;
    m[Y] = 0;
    m[Z] = 0;

    const arma::mat & matR0 = m_particles.r0();

    const double x_a = matR0(a, X);
    const double y_a = matR0(a, Y);
    const double z_a = matR0(a, Z);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_a);

    double k_one[3];
    double k_two[3];


    for(int i=0; i<3; i++) {
        k_one[X] = 0;
        k_one[Y] = 0;
        k_one[Z] = 0;

        k_two[X] = 0;
        k_two[Y] = 0;
        k_two[Z] = 0;

        for(auto &con:PDconnections) {
            int id_b = con.first;
            int b = m_idToCol.at(id_b);

            dR0[X] = x_a - matR0(b, X);
            dR0[Y] = y_a - matR0(b, Y);
            dR0[Z] = z_a - matR0(b, Z);

            const double c_b = m_data(b, m_indexMicromodulus);
            const double vol_b = m_data(b, m_indexVolume);
            const double dr0Len = con.second[m_indexDr0];
            const double volumeScaling = con.second[m_indexVolumeScaling];
            const double g_ab = con.second[m_indexForceScaling];
            const double c_ab = 0.5*(c_a + c_b)*g_ab;

            const double dr0Len3 = pow(dr0Len, 3);
            const double coeff = c_ab*volumeScaling*vol_b/dr0Len3;
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

        m[i] = k_two[i];
    }

    double stiffness = 0;

    for(int d=0;d<3; d++) {
        if(m[d]>stiffness) {
            stiffness = m[d];
        }
    }

    return 0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void PD_PMB::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);


    // Setting the initial max stretch between two particles
    const ivec & colToId = m_particles.colToId();
    for(unsigned int i=0;i<m_particles.nParticles();i++) {
        const int pId = colToId(i);
        const double s0_i = m_data(i, m_indexS0);

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int col_j = m_idToCol.at(id_j);
            const double s0_j = m_data(col_j, m_indexS0);
            const double s0 = 0.5*(s0_i + s0_j);
            con.second[m_indexS00] = s0;
            m_data(i, m_indexS_new) = s0;
        }
    }

    // Setting the micromodulus
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
        nu = 1.;
        k = E;
        c = 12*E/(h*h*pow(delta, 2));
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

    const double dimScaling = 2.*pow(dim, 2.);

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++) {
        const int pId = colToId(i);
        double dRvolume = 0;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int col_j = m_idToCol.at(id_j);
            const double dr0Len = con.second[m_indexDr0];

            dRvolume += dr0Len*m_data(col_j, m_indexVolume);
//            double volumeScaling = con.second[m_indexVolumeScaling];
            //            dRvolume += dr0Len*data(col_j, indexVolume)*volumeScaling;
        }
        m_data(i, m_indexMicromodulus) *= dimScaling/dRvolume;
    }
}
//------------------------------------------------------------------------------
}

