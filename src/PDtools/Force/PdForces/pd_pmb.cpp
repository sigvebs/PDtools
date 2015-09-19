#include "pd_pmb.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_PMB::PD_PMB(PD_Particles &particles, double lc, double delta, double alpha):
    Force(particles),
    m_r(m_particles.r()),
    m_r0(m_particles.r0()),
    m_F(m_particles.F()),
    m_data(m_particles.data()),
    m_pIds(m_particles.pIds())
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_indexS0 = m_particles.getParamId("s0");
    m_indexS00 = m_particles.registerPdParameter("s00");
    m_indexConnected = m_particles.getPdParamId("connected");
    m_indexS_new = m_particles.registerParameter("s_new");

    // Setting the initial max stretch between two particles
    for(int i=0;i<m_particles.nParticles();i++)
    {
        int pId = i;
        int col_i = i;
        double s0_i = m_data(col_i, m_indexS0);

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = m_pIds[id_j];
            double s0_j = m_data(col_j, m_indexS0);
            double s0 = 0.5*(s0_i + s0_j);
            con.second[m_indexS00] = s0;
            m_data(i, m_indexS_new) = s0;
        }
    }


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

    for(int i=0;i<nParticles; i++)
    {
        f[i] = m_F.colptr(i);
        x[i] = m_r.colptr(i);
        r0[i] = m_r0.colptr(i);
    }
}
//------------------------------------------------------------------------------
PD_PMB::~PD_PMB()
{
}
//------------------------------------------------------------------------------
void PD_PMB::calculateForces(const pair<int, int> &idCol)
{
    int pId = idCol.first;
    int i = idCol.second;
    double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
//    vector<pair<int, vector<double>> *> removeParticles;

    int jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq, r, dr, rk, fbond;
    double delta, stretch, s0_new;

    jnum = PDconnections.size();
    double lc = m_lc;
    delta = m_delta;
    double half_lc = 0.5*lc;
    double vfrac_scale = 1.0;

    bool first;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    m_data(i, m_indexS0) = m_data(i, m_indexS_new);
    double s0_i = m_data(i, m_indexS0);
    s0_new = numeric_limits<double>::max();
    first = true;

    for (int jj = 0; jj < jnum; jj++)
    {
        auto &con = PDconnections[jj];
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        int id_j = con.first;
        int j = m_pIds[id_j];
        double dr0 = con.second[m_indexDr0];

        // compute force density, add to PD equation of motion
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);
        dr = r - dr0;

        // avoid roundoff errors
        if (fabs(dr) < 2.2204e-016)
        {
            dr = 0.0;
        }

        // scale vfrac[j] if particle j near the horizon
        double c_j = m_data(j, m_indexMicromodulus);
        double c_ij = 0.5*(c_i + c_j);

        if ((fabs(dr0 - delta)) <= half_lc)
        {
            vfrac_scale = (-1.0/(2*half_lc))*(dr0)
                    + (1.0 + ((delta - half_lc)/(2*half_lc) ) );
        }
        else
        {
            vfrac_scale = 1.0;
        }

        double vfrac_j = m_data(j, m_indexVolume);

        stretch = dr / dr0;
        rk = (c_ij * vfrac_j) * vfrac_scale * stretch;

        if (r > 0.0)
            fbond = -(rk/r);
        else
            fbond = 0.0;

        f[i][0] += delx*fbond;
        f[i][1] += dely*fbond;
        f[i][2] += delz*fbond;

        double s0_j = m_data(j, m_indexS0);

        if (stretch > min(s0_i, s0_j)){
            con.second[m_indexConnected] = 0;
        }
//            removeParticles.push_back(&con);

        // update s0 for next timestep
        double s00 = con.second[m_indexS00];
        s0_new = s00;
        /*
        if (first)
           s0_new = s00 - (m_alpha * stretch);
        else
           s0_new = max(s0_new, s00 - m_alpha * stretch);
        */
        first = false;
      }

    m_data(i, m_indexS_new) = s0_new;
//    for(auto &removeParticle:removeParticles)
//    {
//        PDconnections.erase( remove(begin(PDconnections), end(PDconnections), *removeParticle),
//                             end(PDconnections) );
//    }

    // store new s0
    //    for (i = 0; i < nlocal; i++)
    //          s0[i] = s0_new[i];


    /*
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
        double g_ij = con.second[m_indexForceScaling];
        double c_ij = 0.5*(c_i + c_j)*g_ij;

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
    */
}
//------------------------------------------------------------------------------
double PD_PMB::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
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
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        int id_j = con.first;
        int col_j = m_pIds[id_j];

        double vol_j = m_data(col_j, m_indexVolume);
        double dr0Len         = con.second[m_indexDr0];
        double volumeScaling   = con.second[m_indexVolumeScaling];
        double c_j = m_data(col_j, m_indexMicromodulus);
        double g_ij = con.second[m_indexForceScaling];

        double c_ij = 0.5*(c_i + c_j)*g_ij;


        dr_ij[X] = m_r(X, col_j) - x_i;
        dr_ij[Y] = m_r(Y, col_j) - y_i;
        dr_ij[Z] = m_r(Z, col_j) - z_i;

        double drSquared = dr_ij[X]*dr_ij[X] + dr_ij[Y]*dr_ij[Y] + dr_ij[Z]*dr_ij[Z];
        double drLen = sqrt(drSquared);
        double ds = drLen - dr0Len;

        // To avoid roundoff errors
        if (fabs(ds) < THRESHOLD)
            ds = 0.0;

        energy += c_ij*(ds*ds)/dr0Len*vol_j*volumeScaling;
    }

    return 0.25*energy;
}
//------------------------------------------------------------------------------
void PD_PMB::calculatePotentialEnergy(const std::pair<int, int> &idCol, int indexPotential)
{
    int col_i = idCol.second;
    double vol_i = m_data(col_i, m_indexVolume);
    m_data(col_i, indexPotential) += calculatePotentialEnergyDensity(idCol)*vol_i;
}
//------------------------------------------------------------------------------
void PD_PMB::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    int pId = idCol.first;
    int i = idCol.second;
    double c_i = m_data(i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
    vector<pair<int, vector<double>> *> removeParticles;

    // loop over my particles and their partners
    // partner list contains all bond partners, so I-J appears twice
    // if bond already broken, skip this partner
    // first = true if this is first neighbor of particle i
    int jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq, r, dr, rk, fbond;
    double delta, stretch, s0_new;

    jnum = PDconnections.size();
    double lc = m_lc;
    delta = m_delta;
    double half_lc = 0.5*lc;
    double vfrac_scale = 1.0;

    bool first;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    double f_tmp[3];

    m_data(i, m_indexS0) = m_data(i, m_indexS_new);
    s0_new = numeric_limits<double>::max();
    first = true;

    for (int jj = 0; jj < jnum; jj++) {
        auto &con = PDconnections[jj];
        if(con.second[m_indexConnected] <= 0.5)
            continue;
        int id_j = con.first;
        int j = m_pIds[id_j];
        double dr0 = con.second[m_indexDr0];

        // compute force density, add to PD equation of motion
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);
        dr = r - dr0;

        // avoid roundoff errors
        if (fabs(dr) < 2.2204e-016)
        {
            dr = 0.0;
        }

        // scale vfrac[j] if particle j near the horizon
        double c_j = m_data(j, m_indexMicromodulus);
        double c_ij = 0.5*(c_i + c_j);

        if ((fabs(dr0 - delta)) <= half_lc)
        {
            vfrac_scale = (-1.0/(2*half_lc))*(dr0)
                    + (1.0 + ((delta - half_lc)/(2*half_lc) ) );
        }
        else
        {
            vfrac_scale = 1.0;
        }

        double vfrac_j = m_data(j, m_indexVolume);

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
        m_data(i, indexStress[2]) += 0.5*f_tmp[Z]*delz;
        m_data(i, indexStress[3]) += 0.5*f_tmp[X]*dely;
        m_data(i, indexStress[4]) += 0.5*f_tmp[X]*delz;
        m_data(i, indexStress[5]) += 0.5*f_tmp[Z]*delz;
      }
}
//------------------------------------------------------------------------------
double PD_PMB::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
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
            double g_ab = con.second[m_indexForceScaling];
            double c_ab = 0.5*(c_a + c_b)*g_ab;

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
}
//------------------------------------------------------------------------------
void PD_PMB::initialize(double E, double nu, double delta, int dim, double h)
{
    Force::initialize(E, nu, delta, dim, h);
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
        nu = 1.;
        k = E;
        c = 12*E/(h*h*pow(delta, 2));
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

