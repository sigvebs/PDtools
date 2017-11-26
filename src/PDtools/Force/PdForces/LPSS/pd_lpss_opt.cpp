#include "pd_lpss_opt.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_LPSS_opt::PD_LPSS_opt(PD_Particles &particles, bool planeStress):
    PD_LPSS(particles, planeStress)
{
//    m_hasUpdateState = true;
//    m_hasStepOneModifier = true;
    //----------------------------------
    // Registering the rotation matrix and the shape matrix
    if(m_dim == 2) {
        m_iRn[0] = m_particles.registerParameter("Rn_00", 1);
        m_iRn[1] = m_particles.registerParameter("Rn_10", 0);
        m_iRn[2] = m_particles.registerParameter("Rn_01", 0);
        m_iRn[3] = m_particles.registerParameter("Rn_11", 1);
    } else if(m_dim == 3) {
        m_iRn[0] = m_particles.registerParameter("Rn_00", 1);
        m_iRn[1] = m_particles.registerParameter("Rn_10", 0);
        m_iRn[2] = m_particles.registerParameter("Rn_20", 0);
        m_iRn[3] = m_particles.registerParameter("Rn_01", 0);
        m_iRn[4] = m_particles.registerParameter("Rn_11", 1);
        m_iRn[5] = m_particles.registerParameter("Rn_21", 0);
        m_iRn[6] = m_particles.registerParameter("Rn_02", 0);
        m_iRn[7] = m_particles.registerParameter("Rn_12", 0);
        m_iRn[8] = m_particles.registerParameter("Rn_22", 1);
    }
}
//------------------------------------------------------------------------------
void PD_LPSS_opt::calculateForces(const int id, const int i)
{
    const double theta_i = m_theta[i];
    const double m_i = m_mass[i];

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);
    const int nConnections = PDconnections.size();

    double dr_ij[m_dim];
    double dr0_ij[m_dim];
    double drr_ij[m_dim];
    double drr_ji[m_dim];
    double dud_ij[m_dim];
    double dud_ji[m_dim];

    int nConnected = 0;
    mat F = zeros(m_dim, m_dim);
    mat K = zeros(m_dim, m_dim);

    double theta_new = 0;

    double R1[5*m_dim - 6]; // 4 for d=2 and 9 for d=3
    double R2[5*m_dim - 6];

    if(m_dim == 3) {
        R1[0] = m_data(i, m_iR[0]); // 00
        R1[1] = m_data(i, m_iR[1]); // 10
        R1[2] = m_data(i, m_iR[2]); // 20
        R1[3] = m_data(i, m_iR[3]); // 01
        R1[4] = m_data(i, m_iR[4]); // 11
        R1[5] = m_data(i, m_iR[5]); // 21
        R1[6] = m_data(i, m_iR[6]); // 02
        R1[7] = m_data(i, m_iR[7]); // 12
        R1[8] = m_data(i, m_iR[8]); // 22

        K(0,0) = m_data(i, m_iK[0]);
        K(1,1) = m_data(i, m_iK[1]);
        K(0,1) = m_data(i, m_iK[2]);
        K(1,0) = K(0,1);
        K(2,2) = m_data(i, m_iK[3]);
        K(0,2) = m_data(i, m_iK[4]);
        K(2,0) = K(0,2);
        K(1,2) = m_data(i, m_iK[5]);
        K(2,1) = K(1,2);

        for(int i=0; i<6; i++) {
            m_data(i, m_iStress[i]) = 0;
        }
        //----------------------------------
        for(int l_j=0; l_j<nConnections; l_j++) {
            auto &con = PDconnections[l_j];
            vector<double> & con_data = con.second;
            if(con_data[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_idToCol.at(id_j);

            R2[0] = m_data(j, m_iR[0]); // 00
            R2[1] = m_data(j, m_iR[1]); // 10
            R2[2] = m_data(j, m_iR[2]); // 20
            R2[3] = m_data(j, m_iR[3]); // 01
            R2[4] = m_data(j, m_iR[4]); // 11
            R2[5] = m_data(j, m_iR[5]); // 21
            R2[6] = m_data(j, m_iR[6]); // 02
            R2[7] = m_data(j, m_iR[7]); // 12
            R2[8] = m_data(j, m_iR[8]); // 22

            const double m_j = m_mass[j];
            const double theta_j = m_theta[j];
            const double vol_j = m_volume[j];
            const double volumeScaling = con_data[m_iVolumeScaling];
            const double vol = vol_j*volumeScaling;
            const double dr0 = con_data[m_iDr0];
            const double w = weightFunction(dr0);

            dr_ij[0] = m_x[j] - m_x[i];
            dr_ij[1] = m_y[j] - m_y[i];
            dr_ij[2] = m_z[j] - m_z[i];
            dr0_ij[0] = m_x0[j] - m_x0[i];
            dr0_ij[1] = m_y0[j] - m_y0[i];
            dr0_ij[2] = m_z0[j] - m_z0[i];

            // drr_ij = R_i*dr0_ij;
            drr_ij[0] = R1[0]*dr0_ij[0] + R1[3]*dr0_ij[1] + R1[6]*dr0_ij[2];
            drr_ij[1] = R1[1]*dr0_ij[0] + R1[4]*dr0_ij[1] + R1[7]*dr0_ij[2];
            drr_ij[2] = R1[2]*dr0_ij[0] + R1[5]*dr0_ij[1] + R1[8]*dr0_ij[2];

            dud_ij[0] = dr_ij[0] - (1. + theta_i/m_dim)*drr_ij[0];
            dud_ij[1] = dr_ij[1] - (1. + theta_i/m_dim)*drr_ij[1];
            dud_ij[2] = dr_ij[2] - (1. + theta_i/m_dim)*drr_ij[2];

            // drr_ji = -R_j*dr0_ij;
            drr_ji[0] = - R2[0]*dr0_ij[0] - R2[3]*dr0_ij[1] - R2[6]*dr0_ij[2];
            drr_ji[1] = - R2[1]*dr0_ij[0] - R2[4]*dr0_ij[1] - R2[7]*dr0_ij[2];
            drr_ji[2] = - R2[2]*dr0_ij[0] - R2[5]*dr0_ij[1] - R2[8]*dr0_ij[2];

            dud_ji[0] = -dr_ij[0] - (1. + theta_j/m_dim)*drr_ji[0];
            dud_ji[1] = -dr_ij[1] - (1. + theta_j/m_dim)*drr_ji[1];
            dud_ji[2] = -dr_ij[2] - (1. + theta_j/m_dim)*drr_ji[2];

            m_Fx[i] += vol*m_dim*w*((m_k*theta_i*drr_ij[0] + 2*m_mu*dud_ij[0])/m_i - (m_k*theta_j*drr_ji[0] + 2*m_mu*dud_ji[0])/m_j);
            m_Fy[i] += vol*m_dim*w*((m_k*theta_i*drr_ij[1] + 2*m_mu*dud_ij[1])/m_i - (m_k*theta_j*drr_ji[1] + 2*m_mu*dud_ji[1])/m_j);
            m_Fz[i] += vol*m_dim*w*((m_k*theta_i*drr_ij[2] + 2*m_mu*dud_ij[2])/m_i - (m_k*theta_j*drr_ji[2] + 2*m_mu*dud_ji[2])/m_j);

            const double ds = sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2])-dr0;
            con_data[m_iStretch] = ds/dr0;

            // Calculating the new deformation gradient
            for(int d=0; d<m_dim; d++) {
                for(int d2=0; d2<m_dim; d2++) {
                    F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
                }
            }

            //theta_new
            double dudr0 = 0;
            for(int d=0; d<m_dim; d++) {
                dudr0 += (dr_ij[d] - drr_ij[d])*drr_ij[d];
            }
            theta_new += w*dudr0*vol;

            nConnected++;
        }
    }
    if (m_dim == 2) { // dim2
        R1[0] = m_data(i, m_iR[0]); // 00
        R1[1] = m_data(i, m_iR[1]); // 10
        R1[2] = m_data(i, m_iR[2]); // 01
        R1[3] = m_data(i, m_iR[3]); // 11

        K(0,0) = m_data(i, m_iK[0]);
        K(1,1) = m_data(i, m_iK[1]);
        K(0,1) = m_data(i, m_iK[2]);
        K(1,0) = K(0,1);

        for(int i=0; i<3; i++)
            m_data(i, m_iStress[i]) = 0;

        //----------------------------------
        for(int l_j=0; l_j<nConnections; l_j++) {
            auto &con = PDconnections[l_j];
            vector<double> & con_data = con.second;
            if(con_data[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_idToCol.at(id_j);

            R2[0] = m_data(j, m_iR[0]); // 00
            R2[1] = m_data(j, m_iR[1]); // 10
            R2[2] = m_data(j, m_iR[2]); // 01
            R2[3] = m_data(j, m_iR[3]); // 11

            const double m_j = m_mass[j];
            const double theta_j = m_theta[j];
            const double vol_j = m_volume[j];
            const double volumeScaling = con_data[m_iVolumeScaling];
            const double vol = vol_j*volumeScaling;
            const double dr0 = con_data[m_iDr0];
            const double w = weightFunction(dr0);

            dr_ij[0] = m_x[j] - m_x[i];
            dr_ij[1] = m_y[j] - m_y[i];
            dr0_ij[0] = m_x0[j] - m_x0[i];
            dr0_ij[1] = m_y0[j] - m_y0[i];

            // drr_ij = R_i*dr0_ij;
            drr_ij[0] = R1[0]*dr0_ij[0] + R1[2]*dr0_ij[1];
            drr_ij[1] = R1[1]*dr0_ij[0] + R1[3]*dr0_ij[1];
            dud_ij[0] = dr_ij[0] - (1. + theta_i/m_dim)*drr_ij[0];
            dud_ij[1] = dr_ij[1] - (1. + theta_i/m_dim)*drr_ij[1];
            // dud_ij = dr_ij - (1. + theta_i/m_dim)*drr_ij;

            // drr_ji = -R_j*dr0_ij;
            drr_ji[0] = - R2[0]*dr0_ij[0] - R2[2]*dr0_ij[1];
            drr_ji[1] = - R2[1]*dr0_ij[0] - R2[3]*dr0_ij[1];
            dud_ji[0] = -dr_ij[0] - (1. + theta_j/m_dim)*drr_ji[0];
            dud_ji[1] = -dr_ij[1] - (1. + theta_j/m_dim)*drr_ji[1];
            // dud_ji = -dr_ij - (1. + theta_j/m_dim)*drr_ji;

            m_Fx[i] += vol*m_dim*w*((m_k*theta_i*drr_ij[0] + 2*m_mu*dud_ij[0])/m_i - (m_k*theta_j*drr_ji[0] + 2*m_mu*dud_ji[0])/m_j);
            m_Fy[i] += vol*m_dim*w*((m_k*theta_i*drr_ij[1] + 2*m_mu*dud_ij[1])/m_i - (m_k*theta_j*drr_ji[1] + 2*m_mu*dud_ji[1])/m_j);

            // Setting the stretch
            const double ds = sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1])-dr0;
            con_data[m_iStretch] = ds/dr0;

            // Calculating the new deformation gradient
            for(int d=0; d<m_dim; d++) {
                for(int d2=0; d2<m_dim; d2++) {
                    F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
                }
            }

//            const double dudr0 = (dr_ij[0] - drr_ij[0])*drr_ij[0] + (dr_ij[1] - drr_ij[1])*drr_ij[1];
//            theta_new += w*dudr0*vol;

            nConnected++;
        }
    }

//    theta_new = m_dim/m_i*theta_new;
//    m_data(i, m_iThetaNew) = theta_new;
    // Computing the new rotational matrix
    F *= K;

    mat U(m_dim, m_dim);
    vec s(m_dim);
    mat V(m_dim, m_dim);
    svd(U, s, V, F);
    mat R = U*V.t();

    if(m_dim == 2) {
        m_data(i, m_iRn[0]) = R(0,0); // 00
        m_data(i, m_iRn[1]) = R(1,0); // 10
        m_data(i, m_iRn[2]) = R(0,1); // 01
        m_data(i, m_iRn[3]) = R(1,1); // 11
    }
    if(m_dim == 3) {
        m_data(i, m_iRn[0]) = R(0,0); // 00
        m_data(i, m_iRn[1]) = R(1,0); // 10
        m_data(i, m_iRn[2]) = R(2,0); // 20
        m_data(i, m_iRn[3]) = R(0,1); // 01
        m_data(i, m_iRn[4]) = R(1,1); // 11
        m_data(i, m_iRn[5]) = R(2,1); // 21
        m_data(i, m_iRn[6]) = R(0,2); // 02
        m_data(i, m_iRn[7]) = R(1,2); // 12
        m_data(i, m_iRn[8]) = R(2,2); // 22
    }
    m_continueState = false;

    /*
    double theta = 0;
    vec dr0_ij_vec = zeros(m_dim);
    vec drr_ij_vec = zeros(m_dim);
    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con = PDconnections[l_j];
        const vector<double> & con_data = con.second;

        if(con_data[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling_ij = con.second[m_iVolumeScaling];
        const double vol = vol_j*volumeScaling_ij;
        const double w = weightFunction(dr0);

        for(int d=0; d<m_dim; d++) {
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr0_ij_vec[d] = m_r0(j, d) - m_r0(i, d);
        }
        drr_ij_vec = R*dr0_ij_vec;
        double dudr0 = 0;
        for(int d=0; d<m_dim; d++) {
            dudr0 += (dr_ij[d] - drr_ij_vec[d])*drr_ij_vec[d];
        }
        theta += w*dudr0*vol;
    }
    theta *= m_dim/m_i;
    m_data(i, m_iThetaNew) = theta;
    */

}
//------------------------------------------------------------------------------
void PD_LPSS_opt::evaluateStepOne()
{
    const int nParticles = m_particles.nParticles();

    // Checking for broken state
    /*
    for(int i=0; i<nParticles; i++) {
        const int id_i = m_colToId.at(i);
        const int broken  = m_data(i, m_iBrokenNow);
        if(broken) {
            computeMandK(id_i, i);
            m_data(i, m_iBrokenNow) = 0;
        }
    }
    */
#if 1
//    mat F = zeros(m_dim, m_dim);
//    mat K = zeros(m_dim, m_dim);
    mat R = zeros(m_dim, m_dim);
//    mat U = zeros(m_dim, m_dim);
//    vec s = zeros(m_dim);
//    mat V = zeros(m_dim, m_dim);

    vec dr_ij = zeros(m_dim);
    vec dr0_ij = zeros(m_dim);
    vec drr_ij = zeros(m_dim);

    for(int i=0; i<nParticles; i++) {
        const int id = m_colToId.at(i);

        const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);
        const int nConnections = PDconnections.size();

        if(m_dim == 2) {
            m_data(i, m_iR[0]) = m_data(i, m_iRn[0]); // 00
            m_data(i, m_iR[1]) = m_data(i, m_iRn[1]); // 10
            m_data(i, m_iR[2]) = m_data(i, m_iRn[2]); // 01
            m_data(i, m_iR[3]) = m_data(i, m_iRn[3]); // 11
//            K(0,0) = m_data(i, m_iK[0]);
//            K(1,1) = m_data(i, m_iK[1]);
//            K(0,1) = m_data(i, m_iK[2]);
//            K(1,0) = K(0,1);
            R[0] = m_data(i, m_iR[0]); // 00
            R[1] = m_data(i, m_iR[1]); // 10
            R[2] = m_data(i, m_iR[2]); // 01
            R[3] = m_data(i, m_iR[3]); // 11
        }
        if(m_dim == 3) {
            m_data(i, m_iR[0]) = m_data(i, m_iRn[0]); // 00
            m_data(i, m_iR[1]) = m_data(i, m_iRn[1]); // 10
            m_data(i, m_iR[2]) = m_data(i, m_iRn[2]); // 20
            m_data(i, m_iR[3]) = m_data(i, m_iRn[3]); // 01
            m_data(i, m_iR[4]) = m_data(i, m_iRn[4]); // 11
            m_data(i, m_iR[5]) = m_data(i, m_iRn[5]); // 21
            m_data(i, m_iR[6]) = m_data(i, m_iRn[6]); // 02
            m_data(i, m_iR[7]) = m_data(i, m_iRn[7]); // 12
            m_data(i, m_iR[8]) = m_data(i, m_iRn[8]); // 22
//            K(0,0) = m_data(i, m_iK[0]);
//            K(1,1) = m_data(i, m_iK[1]);
//            K(0,1) = m_data(i, m_iK[2]);
//            K(1,0) = K(0,1);
//            K(2,2) = m_data(i, m_iK[3]);
//            K(0,2) = m_data(i, m_iK[4]);
//            K(2,0) = K(0,2);
//            K(1,2) = m_data(i, m_iK[5]);
//            K(2,1) = K(1,2);
        }

        //----------------------------------
        // Calculating the deformation gradient, F,
        // and the rotation matrix, R.
        //----------------------------------
        /*
        F.zeros();
        for(int l_j=0; l_j<nConnections; l_j++) {
            auto &con = PDconnections[l_j];
            const vector<double> & con_data = con.second;

            if(con_data[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_idToCol.at(id_j);

            const double vol_j = m_data(j, m_iVolume);
            const double dr0 = con.second[m_iDr0];
            const double volumeScaling_ij = con.second[m_iVolumeScaling];
            const double vol = vol_j*volumeScaling_ij;
            const double w = weightFunction(dr0);

            // Computing the deformation matrix
            for(int d=0; d<m_dim; d++) {
                dr_ij[d] = m_r(j, d) - m_r(i, d);
                dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
            }
            for(int d=0; d<m_dim; d++) {
                for(int d2=0; d2<m_dim; d2++) {
                    F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
                }
            }
        }
        F *= K;

        svd(U, s, V, F);
        R = U*V.t();

        if(m_dim >= 2) {
            m_data(i, m_iR[0]) = R(0,0);
            m_data(i, m_iR[1]) = R(1,0);
            m_data(i, m_iR[2]) = R(0,1);
            m_data(i, m_iR[3]) = R(1,1);
        }
        if(m_dim == 3) {
            m_data(i, m_iR[0]) = R(0,0); // 00
            m_data(i, m_iR[1]) = R(1,0); // 10
            m_data(i, m_iR[2]) = R(2,0); // 20
            m_data(i, m_iR[3]) = R(0,1); // 01
            m_data(i, m_iR[4]) = R(1,1); // 11
            m_data(i, m_iR[5]) = R(2,1); // 21
            m_data(i, m_iR[6]) = R(0,2); // 02
            m_data(i, m_iR[7]) = R(1,2); // 12
            m_data(i, m_iR[8]) = R(2,2); // 22
        }
        */
        //----------------------------------
        // Calculating the dilation, theta
        //----------------------------------
        const double m_i = m_mass[i];
        double theta = 0;

        for(int l_j=0; l_j<nConnections; l_j++) {
            auto &con = PDconnections[l_j];
            const vector<double> & con_data = con.second;

            if(con_data[m_iConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = m_idToCol.at(id_j);

            const double vol_j = m_data(j, m_iVolume);
            const double dr0 = con.second[m_iDr0];
            const double volumeScaling_ij = con.second[m_iVolumeScaling];
            const double vol = vol_j*volumeScaling_ij;
            const double w = weightFunction(dr0);

            for(int d=0; d<m_dim; d++) {
                dr_ij[d] = m_r(j, d) - m_r(i, d);
                dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
            }
            drr_ij = R*dr0_ij;
            double dudr0 = 0;
            for(int d=0; d<m_dim; d++) {
                dudr0 += (dr_ij[d] - drr_ij[d])*drr_ij[d];
            }
            theta += w*dudr0*vol;
        }
        theta *= m_dim/m_i;
        m_theta[i] = theta;
    }
#endif
}
//------------------------------------------------------------------------------
void PD_LPSS_opt::evaluateStepOne(const int id_i, const int i)
{
    (void) id_i;
//    m_theta[i] =  m_data(i, m_iThetaNew);
    /*
    */
//    m_theta[i] =  m_data(i, m_iThetaNew);

}
//------------------------------------------------------------------------------
}
