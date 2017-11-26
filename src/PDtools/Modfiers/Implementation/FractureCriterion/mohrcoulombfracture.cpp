#include "mohrcoulombfracture.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
MohrCoulombFracture::MohrCoulombFracture(double mu, double C, double T):
    m_S0(C), m_T(T)
{
    m_neededProperties= {pair<string, int>("stress", 1)};
    m_phi = mu*M_PI/180.;
    m_d = tan(m_phi);

    m_cos_theta = cos(M_PI/2. + m_phi);
    m_sin_theta = sin(M_PI/2. + m_phi);

    m_S0 = C;
    m_C0 = 2.*m_S0*m_cos_theta/(m_sin_theta - 1.);
    m_k = (m_sin_theta + 1.)/(m_sin_theta - 1.);

    m_hasStepOne = true;
}
//------------------------------------------------------------------------------
void MohrCoulombFracture::initialize()
{

}
//------------------------------------------------------------------------------
void MohrCoulombFracture::registerParticleParameters()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_particles->registerParameter("damage");
    m_idToCol = &m_particles->idToCol();
    m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);

    switch(m_dim)
    {
    case 1:
        m_ghostParameters = {"s_xx"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        break;
    case 2:
        m_ghostParameters = {"s_xx", "s_yy", "s_xy"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        break;
    case 3:
        m_ghostParameters = {"s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStress[3] = m_particles->registerParameter("s_zz");
        m_indexStress[4] = m_particles->registerParameter("s_xz");
        m_indexStress[5] = m_particles->registerParameter("s_yz");
        break;
    }
}
//------------------------------------------------------------------------------
void MohrCoulombFracture::evaluateStepOne(const int id_i, const int i)
{
    // First calculating the total stress on an material point. The stress is averaged
    // all its bonds, and the stress state on a bond is the mean of the
    // stress at each material point.

    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;
    mat & data = *m_data;

#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
    arma::vec eigval(m_dim);
#endif
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    if(m_dim == 2) {
        for(auto &con:PDconnections) {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];

            if((*m_data)(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const double sx = 0.5*(data(i, m_indexStress[0]) + data(j, m_indexStress[0]));
            const double sy = 0.5*(data(i, m_indexStress[1]) + data(j, m_indexStress[1]));
            const double sxy = 0.5*(data(i, m_indexStress[2]) + data(j, m_indexStress[2]));

            const double first = 0.5*(sx + sy);
            const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

            const double p_1 = first + second; // max
            const double p_2 = first - second; // min

            //            const int criticalShear = p_2 <= m_C0 - m_k*p_1;
            const int criticalShear = 0;
            const int criticalTensile = p_1 >= m_T;

            if(criticalShear) {
                data(i, m_indexBrokenNow) = 1;
                con.second[m_indexConnected] = 0;
            } else if(criticalTensile) {
                data(i, m_indexBrokenNow) = 1;
                con.second[m_indexConnected] = 0;
            }

            //            double p_2 = first + second; // max
            //            double p_1 = first - second; // min
//            const double shear = fabs(0.5*(p_1 - p_2)*m_sin_theta);
//            const double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*m_cos_theta;

//            const double criticalShear = fabs(shear) - fabs(m_S0 - m_d*normal);
//            const double criticalTensile = p_2 - m_T;

//            if(criticalShear >= 0  && normal < 0)
//            {
//                data(i, m_indexBrokenNow) = 1;
//                con.second[m_indexConnected] = 0;
//            }
//            else if(criticalTensile >= 0)
//            {
//                data(i, m_indexBrokenNow) = 1;
//                con.second[m_indexConnected] = 0;
//            }
        }
    }
    else if(m_dim == 3)
    {
        /*
        arma::mat S_i(m_dim, m_dim);
        arma::mat S(m_dim, m_dim);
        S_i(0, 0) = (*m_data)(i, m_indexStress[0]);
        S_i(1, 1) = (*m_data)(i, m_indexStress[1]);
        S_i(0, 1) = (*m_data)(i, m_indexStress[2]);
        S_i(1, 0) = S_i(0, 1);
        S_i(2, 2) = (*m_data)(i, m_indexStress[3]);
        S_i(0, 2) = (*m_data)(i, m_indexStress[4]);
        S_i(2, 0) = S_i(0, 2);
        S_i(1, 2) = (*m_data)(i, m_indexStress[5]);
        S_i(2, 1) = S_i(1, 2);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];

            if((*m_data)(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(S_i(0, 0) + (*m_data)(j, m_indexStress[0]));
            S(1, 1) = 0.5*(S_i(1, 1) + (*m_data)(j, m_indexStress[1]));
            S(0, 1) = 0.5*(S_i(0, 1) + (*m_data)(j, m_indexStress[2]));
            S(1, 0) = S(0, 1);
            S(2, 2) = 0.5*(S_i(2, 2) + (*m_data)(j, m_indexStress[3]));
            S(0, 2) = 0.5*(S_i(0, 2) + (*m_data)(j, m_indexStress[4]));
            S(2, 0) = S(0, 2);
            S(1, 2) = 0.5*(S_i(1, 2) + (*m_data)(j, m_indexStress[5]));
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
            const double s2 = I1/3. + core*cos(phi + 2.*M_PI/3.);
            const double s3 = I1/3. + core*cos(phi + 4.*M_PI/3.);

            double p_1 = s1;
            double p_2 = s2;

            if(s2>p_1)
            {
                p_1 = s2;
                p_2 = s1;
            }
            if(s3>p_1)
            {
                p_1 = s3;
            }
            else if(p_2 < s3){
                p_2 = s3;
            }
#endif
            if(m_d*p_1 - p_2 - m_C > 0)
            {
                con.second[m_indexConnected] = 0;
            }
            else
            if(p_1 > m_T)
            {
                con.second[m_indexConnected] = 0;
            }
        }*/
    }
}
//------------------------------------------------------------------------------
}

