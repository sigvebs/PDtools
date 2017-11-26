#include "adrmohrcoulombfracture.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ADRmohrCoulombFracture::ADRmohrCoulombFracture(double mu, double C, double T):
    m_C(C), m_T(T)
{
    m_neededProperties= {pair<string, int>("stress",1)};
    m_phi = mu*M_PI/180.;
    m_d = tan(m_phi);

    m_cos_theta = cos(M_PI/2. + m_phi);
    m_sin_theta = sin(M_PI/2. + m_phi);
    m_tan_theta = tan(0.5*(M_PI + m_phi));
    m_tan2_theta = m_tan_theta*m_tan_theta;

    m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::initialize()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_particles->registerParameter("damage");
    m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
    m_idToCol = &m_particles->idToCol();

    m_state = false;
    m_maxPId = pair<int, int>(-1, -1);
    m_maxStress = std::numeric_limits<double>::min();

    switch(m_dim)
    {
    case 1:
        m_nStressElements = 1;
        m_ghostParameters = {"s_xx"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        break;
    case 2:
        m_nStressElements = 3;
        m_ghostParameters = {"s_xx", "s_yy", "s_xy"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        break;
    case 3:
        m_nStressElements = 6;
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
void ADRmohrCoulombFracture::evaluateStepTwo(const int id_i, const int i)
{
    // First calculating the total stress on an material point. The stress is averaged
    // all its bonds, and the stress state on a bond is the mean of the
    // stress at each material point.
    if((*m_data)(i, m_indexUnbreakable) >= 1)
        return;
    mat &data = *m_data;
    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

//    arma::vec eigval(m_dim);
//    arma::mat S_i(m_dim, m_dim);
//    arma::mat S(m_dim, m_dim);

    if(m_dim == 2) {
//        S_i(0, 0) = data(i, m_indexStress[0]);
//        S_i(1, 1) = data(i, m_indexStress[1]);
//        S_i(0, 1) = data(i, m_indexStress[2]);
//        S_i(1, 0) = S_i(0, 1);
        int counter = 0;
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

            const double shear_max = 0.5*(p_1 - p_2);
            const double shear = shear_max*m_cos_theta;
            const double normal = 0.5*(p_1 + p_2) + shear_max*m_sin_theta;

            const double criticalShear = shear + m_d*normal - m_C;
            const double criticalTensile = p_1 - m_T;

            //                if(criticalShear >= 0  && normal < 0)
            if(criticalShear >= 0) {
                con.second[m_indexConnected] = 0;
                data(i, m_indexBrokenNow) = 1;
                m_maxPId = pair<int, int>(id_i, counter);
            } else if(criticalTensile >= 0) {
                con.second[m_indexConnected] = 0;
                data(i, m_indexBrokenNow) = 1;
                m_maxPId = pair<int, int>(id_i, counter);
            }

//            const int id_j = con.first;
//            const int j = (*m_idToCol)[id_j];

//            if(data(j, m_indexUnbreakable) >= 1)
//                continue;

//            if(con.second[m_indexConnected] <= 0.5)
//                continue;

//            S(0, 0) = 0.5*(S_i(0, 0) + data(j, m_indexStress[0]));
//            S(1, 1) = 0.5*(S_i(1, 1) + data(j, m_indexStress[1]));
//            S(0, 1) = 0.5*(S_i(0, 1) + data(j, m_indexStress[2]));
//            S(1, 0) = S(0, 1);
//#if 0
//            arma::eig_sym(eigval, S);
//            const double p_1 = eigval(1);
//            const double p_2 = eigval(0);
//#else
//            const double first = 0.5*(S(0, 0) + S(1, 1));
//            const double second = sqrt(0.25*(S(0, 0) - S(1, 1))*(S(0, 0) - S(1, 1)) + S(0, 1)*S(0, 1));
//            const double s1 = first + second;
//            const double s2 = first - second;
//            const double p_1 = max(s1, s2);
//            const double p_2 = min(s1, s2);
//#endif

//            if(p_1 > m_T)
//            {
//                con.second[m_indexConnected] = 0;
//                m_maxPId = pair<int, int>(id_i, counter);
//            }

//            if(m_d*p_1 - p_2 - m_C > 0)
//            {
//                con.second[m_indexConnected] = 0;
//                m_maxPId = pair<int, int>(id_i, counter);
////                cout << "shearing/compression" << endl;
//            }
//            else
//            if(p_1 > m_T)
//            {
//                con.second[m_indexConnected] = 0;
////                cout << "tension" << endl;
//                if(p_1 > m_maxStress)
//                {
//                    m_maxPId = pair<int, int>(id_i, counter);
//                    m_maxStress = p_1;
//                    con.second[m_indexConnected] = 0;
//                }
//            }
            counter++;
        }
    }
    else if(m_dim == 3) {
        /*
        S_i(0, 0) = data(i, m_indexStress[0]);
        S_i(1, 1) = data(i, m_indexStress[1]);
        S_i(0, 1) = data(i, m_indexStress[2]);
        S_i(1, 0) = S_i(0, 1);
        S_i(2, 2) = data(i, m_indexStress[3]);
        S_i(0, 2) = data(i, m_indexStress[4]);
        S_i(2, 0) = S_i(0, 2);
        S_i(1, 2) = data(i, m_indexStress[5]);
        S_i(2, 1) = S_i(1, 2);

        int counter = 0;
        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_idToCol)[id_j];

            if(data(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(S_i(0, 0) + data(j, m_indexStress[0]));
            S(1, 1) = 0.5*(S_i(1, 1) + data(j, m_indexStress[1]));
            S(0, 1) = 0.5*(S_i(0, 1) + data(j, m_indexStress[2]));
            S(1, 0) = S(0, 1);
            S(2, 2) = 0.5*(S_i(2, 2) + data(j, m_indexStress[3]));
            S(0, 2) = 0.5*(S_i(0, 2) + data(j, m_indexStress[4]));
            S(2, 0) = S(0, 2);
            S(1, 2) = 0.5*(S_i(1, 2) + data(j, m_indexStress[5]));
            S(2, 1) = S(1, 2);

#if 0
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
                m_maxPId = pair<int, int>(id_i, counter);
            }
            else if(p_1 > m_T)
            {
                con.second[m_indexConnected] = 0;
                m_maxPId = pair<int, int>(id_i, counter);
            }
            counter++;
        }
        */
    }
    //--------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepTwo()
{
    if(m_maxPId.first != -1)
    {
        const int id_i = m_maxPId.first;
        const int remove = m_maxPId.second;
        m_state = true;
//        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
//        PDconnections[remove].second[m_indexConnected] = 0;
    }
    else
    {
        m_state = false;
    }

    m_maxPId = std::pair<int, int>(-1, -1);
    m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
}

