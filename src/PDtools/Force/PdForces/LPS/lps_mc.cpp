#include "lps_mc.h"
#include "PDtools/Particles/pd_particles.h"

#define USE_PRINCIPAL_STRESS 1

namespace PDtools
{
//------------------------------------------------------------------------------
LPS_mc::LPS_mc(PD_Particles &particles, double c, double phi, double C, double T, bool planeStress, bool analyticalM):
    PD_LPS(particles, planeStress, analyticalM),
    m_dampCoeff(c), m_T(T)
{
    particles.setNeedGhostVelocity(true);
    m_particles.setNeedGhostR0(true);
    m_hasStepTwoModifier = true;

    m_phi = phi*M_PI/180.;
    m_d = tan(m_phi);

    m_cos_theta = cos(m_phi);
    m_sin_theta = sin(m_phi);
    m_tan_theta = tan(m_phi);

    m_S0 = C;
    m_C0 = 2.*m_S0*m_cos_theta/(m_sin_theta - 1.);
    m_ks = (m_sin_theta + 1.)/(m_sin_theta - 1.);

    m_smallStrain = true;
//    m_greenStrain = true;
}
//------------------------------------------------------------------------------
void LPS_mc::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    PD_LPS::initialize(E, nu, delta, dim, h, lc);

    _F = zeros(m_dim, m_dim);
    m_strain = zeros(m_dim, m_dim);
    m_K = zeros(m_dim, m_dim);
    m_P = zeros(m_dim, m_dim);

    m_indexConnected = m_particles.getPdParamId("connected");
    m_indexUnbreakable =  m_particles.registerParameter("unbreakable");
    m_indexVolume = m_particles.getParamId("volume");
    m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    m_indexDr0 = m_particles.getPdParamId("dr0");
    m_indexBrokenNow = m_particles.registerParameter("brokenNow", 0);
    m_particles.registerParameter("damage");

    switch(m_dim) {
    case 1:
        m_nStressStrainElements = 1;
        m_ghostParameters.push_back("s_xx");
        m_indexStress[0] = m_particles.registerParameter("s_xx");
        m_indexStrain[0] = m_particles.registerParameter("e_xx");
        m_indexK[0] = m_particles.registerParameter("K_x");
        break;
    case 2:
        m_nStressStrainElements = 3;
//        m_ghostParameters.push_back("s_xx");
//        m_ghostParameters.push_back("s_yy");
//        m_ghostParameters.push_back("s_xy");
//        m_indexK[0] = m_particles.registerParameter("K_x");
//        m_indexK[1] = m_particles.registerParameter("K_y");
//        m_indexK[2] = m_particles.registerParameter("K_xy");
//        m_indexStress[0] = m_particles.registerParameter("s_xx");
//        m_indexStress[1] = m_particles.registerParameter("s_yy");
//        m_indexStress[2] = m_particles.registerParameter("s_xy");
//        m_indexStrain[0] = m_particles.registerParameter("e_xx");
//        m_indexStrain[1] = m_particles.registerParameter("e_yy");
//        m_indexStrain[2] = m_particles.registerParameter("e_xy");

        m_ghostParameters.push_back("s2_xx");
        m_ghostParameters.push_back("s2_yy");
        m_ghostParameters.push_back("s2_xy");
        m_indexK[0] = m_particles.registerParameter("K2_x");
        m_indexK[1] = m_particles.registerParameter("K2_y");
        m_indexK[2] = m_particles.registerParameter("K2_xy");
        m_indexStress[0] = m_particles.registerParameter("s2_xx");
        m_indexStress[1] = m_particles.registerParameter("s2_yy");
        m_indexStress[2] = m_particles.registerParameter("s2_xy");
        m_indexStrain[0] = m_particles.registerParameter("e2_xx");
        m_indexStrain[1] = m_particles.registerParameter("e2_yy");
        m_indexStrain[2] = m_particles.registerParameter("e2_xy");
        break;
    case 3:
        m_nStressStrainElements = 6;
        m_ghostParameters.push_back("s_xx");
        m_ghostParameters.push_back("s_yy");
        m_ghostParameters.push_back("s_xy");
        m_ghostParameters.push_back("s_xy");
        m_ghostParameters.push_back("s_xz");
        m_ghostParameters.push_back("s_yz");
        m_indexK[0] = m_particles.registerParameter("K_x");
        m_indexK[1] = m_particles.registerParameter("K_y");
        m_indexK[2] = m_particles.registerParameter("K_xy");
        m_indexK[3] = m_particles.registerParameter("K_z");
        m_indexK[4] = m_particles.registerParameter("K_xz");
        m_indexK[5] = m_particles.registerParameter("K_yz");
        m_indexStress[0] = m_particles.registerParameter("s_xx");
        m_indexStress[1] = m_particles.registerParameter("s_yy");
        m_indexStress[2] = m_particles.registerParameter("s_xy");
        m_indexStress[3] = m_particles.registerParameter("s_zz");
        m_indexStress[4] = m_particles.registerParameter("s_xz");
        m_indexStress[5] = m_particles.registerParameter("s_yz");
        m_indexStrain[0] = m_particles.registerParameter("e_xx");
        m_indexStrain[1] = m_particles.registerParameter("e_yy");
        m_indexStrain[2] = m_particles.registerParameter("e_xy");
        m_indexStrain[3] = m_particles.registerParameter("e_zz");
        m_indexStrain[4] = m_particles.registerParameter("e_xz");
        m_indexStrain[5] = m_particles.registerParameter("e_yz");
        break;
    }

    // Computing the shape tensor
    const ivec &colToId = m_particles.colToId();
    const int nParticles = m_particles.nParticles();

    for(int i=0; i<nParticles; i++) {
        const int id_i = colToId(i);
        computeK(id_i, i);
    }

    m_lambda = m_E*m_nu/((1+m_nu)*(1. - 2.*m_nu));
    m_mu = 0.5*m_E/(1. + m_nu);
}

//------------------------------------------------------------------------------
void LPS_mc::calculateForces(const int id, const int i)
{
    const double theta_i = m_data(i, m_iTheta);
    const double m_i = m_data(i, m_iMass);

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id);

    const int nConnections = PDconnections.size();
    double dr0_ij[m_dim];
    double dr_ij[m_dim];
    _F.zeros();

    double thetaNew = 0;
    int nConnected = 0;

    //----------------------------------
    // TMP - standard stres calc from
//    m_data(i, m_indexStress[0]) = 0;
//    m_data(i, m_indexStress[1]) = 0;
//    m_data(i, m_indexStress[2]) = 0;
    //----------------------------------
    for(int l_j=0; l_j<nConnections; l_j++) {
        auto &con = PDconnections[l_j];

        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int id_j = con.first;
        const int j = m_idToCol.at(id_j);

        const double m_j = m_data(j, m_iMass);
        const double theta_j = m_data(j, m_iTheta);
        const double vol_j = m_data(j, m_iVolume);
        const double dr0 = con.second[m_iDr0];
        const double volumeScaling = con.second[m_iVolumeScaling];
        const double vol = vol_j*volumeScaling;
        const double w = weightFunction(dr0);

        double dr2 = 0;
        double drdv = 0;

        for(int d=0; d<m_dim; d++) {
            dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
            dr_ij[d] = m_r(j, d) - m_r(i, d);
            dr2 += dr_ij[d]*dr_ij[d];
            drdv += dr_ij[d]*(m_v(j, d) - m_v(i, d));
        }

        const double dr = sqrt(dr2);
        const double dsdt = drdv/(dr);
        const double ds = dr - dr0 + m_dampCoeff*dsdt;
        double bond = m_c*(theta_i/m_i + theta_j/m_j)*dr0;
        bond += m_alpha*(1./m_i + 1./m_j)*ds;
        bond *= w*vol/dr;
        thetaNew += w*dr0*ds*vol;

        for(int d=0; d<m_dim; d++) {
            m_F(i, d) += dr_ij[d]*bond;

            for(int d2=0; d2<m_dim; d2++) {
                _F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
            }
        }

        //----------------------------------
        // TMP - standard stres calc from
//        m_data(i, m_indexStress[0]) += 0.5*dr_ij[0]*dr_ij[0]*bond;
//        m_data(i, m_indexStress[1]) += 0.5*dr_ij[1]*dr_ij[1]*bond;
//        m_data(i, m_indexStress[2]) += 0.5*dr_ij[0]*dr_ij[1]*bond;
        //----------------------------------

        nConnected++;
    }

    if(nConnections <= 5) {
        m_data(i, m_iThetaNew) = 0;
    }
    else
        m_data(i, m_iThetaNew) = m_dim/m_i*thetaNew;

    //----------------------------------
    // TMP - standard stres calc from
    //----------------------------------
    computeStress(id, i, nConnected);
    //----------------------------------
}
//------------------------------------------------------------------------------
void LPS_mc::evaluateStepTwo(int id_i, int i)
{
    if(m_data(i, m_indexUnbreakable) >= 1)
        return;
#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
    arma::vec eigval(m_dim);
#endif
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    const double shearCrit = m_C0 - m_ks*m_T;

    bool broken = false;
    if(m_dim == 2) {
        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int j = m_idToCol[id_j];

            if(m_data(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const double sx = 0.5*(m_data(i, m_indexStress[0]) + m_data(j, m_indexStress[0]));
            const double sy = 0.5*(m_data(i, m_indexStress[1]) + m_data(j, m_indexStress[1]));
            const double sxy = 0.5*(m_data(i, m_indexStress[2]) + m_data(j, m_indexStress[2]));

            const double first = 0.5*(sx + sy);
            const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

            const double p_1 = first + second; // max
            const double p_2 = first - second; // min

#if USE_PRINCIPAL_STRESS
            const int criticalShear = p_2 <= m_C0 - m_ks*p_1;
#else
            const double shear_max = 0.5*(p_1 - p_2);
            const double shear = shear_max*m_cos_theta;
            const double normal = 0.5*(p_1 + p_2) + shear_max*m_sin_theta;
            const double criticalShear = shear > m_S0 - m_d*normal;
//            const double criticalShear = shear + m_d*normal - m_S0;
//            const double criticalTensile = p_1 - m_T;
#endif
            const int MC_valid = p_2 < shearCrit;
            const int criticalTensile = p_1 >= m_T;

            if (MC_valid) {
                if(criticalShear) {
                    m_data(i, m_indexBrokenNow) = 1;
                    con.second[m_indexConnected] = 0;
                    broken = true;
                }
            } else {
                if(criticalTensile) {
                    m_data(i, m_indexBrokenNow) = 1;
                    con.second[m_indexConnected] = 0;
                    broken = true;
                }
            }
        }
    } else if(m_dim == 3) {
        arma::mat S(m_dim, m_dim);

        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int j = m_idToCol[id_j];

            if(m_data(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(m_data(i, m_indexStress[0]) + m_data(j, m_indexStress[0]));
            S(1, 1) = 0.5*(m_data(i, m_indexStress[1]) + m_data(j, m_indexStress[1]));
            S(0, 1) = 0.5*(m_data(i, m_indexStress[2]) + m_data(j, m_indexStress[2]));
            S(1, 0) = S(0, 1);
            S(2, 2) = 0.5*(m_data(i, m_indexStress[3]) + m_data(j, m_indexStress[3]));
            S(0, 2) = 0.5*(m_data(i, m_indexStress[4]) + m_data(j, m_indexStress[4]));
            S(2, 0) = S(0, 2);
            S(1, 2) = 0.5*(m_data(i, m_indexStress[5]) + m_data(j, m_indexStress[5]));
            S(2, 1) = S(1, 2);

#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
            arma::eig_sym(eigval, S);
            const double p_1 = eigval(2);
            const double p_2 = eigval(0);
#else
            const double I1 = S(0, 0) + S(1, 1) + S(2, 2);
            const double I2 = S(0, 0)*S(1, 1) + S(1, 1)*S(2, 2) + S(3, 3)*S(0, 0) - pow(S(0, 1), 2) - pow(S(1, 2), 2) - pow(S(0, 2), 2);
            const double I3 = S(0, 0)*S(1, 1)*S(2, 2) - S(0, 0)*pow(S(1, 2), 2) - S(1, 1)*pow(S(0, 2), 2) - S(2, 2)*pow(S(0, 1), 2) + 2*S(0, 1)*S(1, 2)*S(0, 2);
            const double phi = 1./3.*acos(0.5*(2*pow(I1, 3) -9*I1*I2 + 27*I3)/pow(pow(I1, 2) - 3*I2, 1.5));

            const double core = 2./3.*(sqrt(I1*I1 - 3*I2));
            const double s1 = I1/3. + core*cos(phi);
            const double s2 = I1/3. + core*cos(phi - 2.*M_PI/3.);
            const double s3 = I1/3. + core*cos(phi - 4.*M_PI/3.);

            double p_1 = s1;
            double p_2 = s2;

            if(s2>p_1) {
                p_1 = s2;
                p_2 = s1;
            }
            if(s3>p_1) {
                p_1 = s3;
            }
            else if(p_2 < s3){
                p_2 = s3;
            }
#endif
            const int criticalShear = p_2 <= m_C0 - m_ks*p_1;
            const int MC_valid = p_2 < shearCrit;
            const int criticalTensile = p_1 >= m_T;

            if (MC_valid) {
                if(criticalShear) {
                    m_data(i, m_indexBrokenNow) = 1;
                    con.second[m_indexConnected] = 0;
                    broken = true;
                }
            } else {
                if(criticalTensile) {
                    m_data(i, m_indexBrokenNow) = 1;
                    con.second[m_indexConnected] = 0;
                    broken = true;
                }
            }
        }
    }

    if(broken) {
        updateWeightedVolume(id_i, i);
        computeK(id_i, i);
//        m_data(i, m_indexBrokenNow) = 0;
    }
}
//------------------------------------------------------------------------------
void LPS_mc::computeStress(const int id, const int i, const int nConnected)
{
    //--------------------------------------------------------------------------
    // Stress
    //--------------------------------------------------------------------------
    vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id);

    if(m_dim >= 2) {
        m_K(0, 0) = m_data(i, m_indexK[0]);
        m_K(1, 1) = m_data(i, m_indexK[1]);
        m_K(0, 1) = m_data(i, m_indexK[2]);
        m_K(1, 0) = m_K(0, 1);
    }
    if(m_dim == 3) {
        m_K(2, 2) = m_data(i, m_indexK[3]);
        m_K(0, 2) = m_data(i, m_indexK[4]);
        m_K(1, 2) = m_data(i, m_indexK[5]);
        m_K(2, 0) = m_K(0, 2);
        m_K(2, 1) = m_K(1, 2);
    }

    if(nConnected <= 5) {
        for(int j=0; j<m_nStressStrainElements; j++) {
            m_data(i, m_indexStrain[j]) = 0;
            m_data(i, m_indexStress[j]) = 0;
        }
        return;
    } else {
        _F = _F*m_K; // K = inv(K);
        if(m_smallStrain) {
            m_strain = 0.5*(_F.t() + _F);
            for(int d=0; d<m_dim;d++) {
                m_strain(d, d) -= 1.;
            }
        } else {
            m_strain = 0.5*_F.t()*_F;
            for(int d=0; d<m_dim;d++) {
                m_strain(d, d) -= 0.5;
            }
        }
    }

    m_data(i, m_indexStrain[0]) = m_strain(0,0);
    m_data(i, m_indexStrain[1]) = m_strain(1,1);
    m_data(i, m_indexStrain[2]) = m_strain(0,1);

    // Assuming linear elasticity
    if(m_dim == 2) {
        // Constituent model, linear elastic
        // Computing the second PK stress
        if(m_planeStress) {
            const double a = m_E/(1. - m_nu*m_nu);
            m_P(0,0) = a*(m_strain(0,0) + m_nu*m_strain(1,1));
            m_P(1,1) = a*(m_strain(1,1) + m_nu*m_strain(0,0));
            m_P(0,1) = a*(1 - m_nu)*m_strain(0,1);
            m_P(1,0) = m_P(0,1);
        } else { // Plane strain
            const double a = m_E/((1. + m_nu)*(1. - 2*m_nu));
            m_P(0,0) = a*((1. - m_nu)*m_strain(0,0) + m_nu*m_strain(1,1));
            m_P(1,1) = a*((1. - m_nu)*m_strain(1,1) + m_nu*m_strain(0,0));
            m_P(0,1) = a*0.5*(1. - 2.*m_nu)*m_strain(0,1);
            m_P(1,0) = m_P(0,1);
        }

        // Converting PK stress to Cauchy stress
        if(m_greenStrain) {
            const double detF = 1./det(_F);
            m_P = detF*_F*m_P*_F.t();
        }
        m_data(i, m_indexStress[0]) = m_P(0,0);
        m_data(i, m_indexStress[1]) = m_P(1,1);
        m_data(i, m_indexStress[2]) = m_P(0,1);
    }

    if(m_dim == 3) {
        m_data(i, m_indexStrain[3]) = m_strain(2,2);
        m_data(i, m_indexStrain[4]) = m_strain(0,2);
        m_data(i, m_indexStrain[5]) = m_strain(1,2);

        m_P(0,0) = (m_lambda + 2*m_mu)*m_strain(0,0) + m_lambda*m_strain(1,1) + m_lambda*m_strain(2,2);
        m_P(1,1) = m_lambda*m_strain(0,0) + (m_lambda + 2*m_mu)*m_strain(1,1) + m_lambda*m_strain(2,2);
        m_P(2,2) = m_lambda*m_strain(0,0) + m_lambda*m_strain(1,1) + (m_lambda + 2*m_mu)*m_strain(2,2);
        m_P(0,1) = 2.*m_mu*m_strain(0,1);
        m_P(0,2) = 2.*m_mu*m_strain(0,2);
        m_P(1,2) = 2.*m_mu*m_strain(1,2);

        // Converting PK2 stress to Cauchy stress
        if(m_greenStrain) {
            const double detF = 1./det(_F);
            m_P = detF*_F*m_P*_F.t();
        }

        m_data(i, m_indexStress[0]) = m_P(0,0);
        m_data(i, m_indexStress[1]) = m_P(1,1);
        m_data(i, m_indexStress[2]) = m_P(0,1);
        m_data(i, m_indexStress[3]) = m_P(2,2);
        m_data(i, m_indexStress[4]) = m_P(0,2);
        m_data(i, m_indexStress[5]) = m_P(1,2);
    }
}
//------------------------------------------------------------------------------
void LPS_mc::computeK(int id, int i)
{
    const auto &idToCol = m_particles.idToCol();
    const mat &r0 = m_particles.r0();
    mat &data = m_particles.data();
    mat K = zeros(m_dim, m_dim);
    vector<pair<int, vector<double>>> & PDconnections_i = m_particles.pdConnections(id);
    const int nConnections = PDconnections_i.size();
    double dr0_ij[m_dim];
    double totalVolume = 0;

    int nConnected = 0;
    for(int l_j=0; l_j<nConnections; l_j++) {
        const auto &con_i = PDconnections_i[l_j];
        const int id_j = con_i.first;
        const int j = idToCol.at(id_j);

        const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
        const double vol_j = data(j, m_indexVolume);
        const double vol = vol_j*volumeScaling_ij;
        totalVolume += vol;

        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const double dr0 = con_i.second[m_indexDr0];
        const double w = weightFunction(dr0);

        for(int d=0; d<m_dim; d++) {
            dr0_ij[d] = r0(j, d) - r0(i, d);
        }

        for(int d=0; d<m_dim; d++) {
            for(int d2=0; d2<m_dim; d2++) {
                K(d, d2) += w*dr0_ij[d]*dr0_ij[d2]*vol;
            }
        }

        nConnected++;
    }

    if(nConnected <= 5) {
        double volumeFrac = totalVolume/nConnections*totalVolume;
        volumeFrac = 1./volumeFrac;
        if(m_dim == 3) {
            K(0,0) = volumeFrac;
            K(1,1) = volumeFrac;
            K(0,1) = 0.;
            K(1,0) = 0.;
            K(2,2) = volumeFrac;
            K(1,2) = 0.;
            K(2,1) = 0.;
        } if(m_dim == 2) {
            K(0,0) = volumeFrac;
            K(1,1) = volumeFrac;
            K(0,1) = 0.;
            K(1,0) = 0.;
        }
    } else {
        K = inv_sympd(K);
    }

    if(m_dim >= 2) {
        data(i, m_indexK[0]) = K(0,0);
        data(i, m_indexK[1]) = K(1,1);
        data(i, m_indexK[2]) = K(0,1);
    }
    if(m_dim == 3) {
        data(i, m_indexK[3]) = K(2,2);
        data(i, m_indexK[4]) = K(0,2);
        data(i, m_indexK[5]) = K(1,2);
    }
}
//------------------------------------------------------------------------------
}
