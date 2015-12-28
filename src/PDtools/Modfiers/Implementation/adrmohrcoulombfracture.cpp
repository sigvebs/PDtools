#include "adrmohrcoulombfracture.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ADRmohrCoulombFracture::ADRmohrCoulombFracture(double mu, double C, double T, int dim):
    m_C(C), m_T(T), m_dim(dim), m_d(pow(sqrt(1 + mu*mu) + mu, 2))
{

}
//------------------------------------------------------------------------------
ADRmohrCoulombFracture::~ADRmohrCoulombFracture()
{

}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::initialize()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_pIds = &m_particles->pIds();

    if(m_particles->hasParameter("s_xx"))
        m_indexStress[0] = m_particles->getParamId("s_xx");
    else
        m_indexStress[0] = m_particles->registerParameter("s_xx");

    if(m_particles->hasParameter("s_yy"))
        m_indexStress[1] = m_particles->getParamId("s_yy");
    else
        m_indexStress[1] = m_particles->registerParameter("s_yy");

    if(m_particles->hasParameter("s_zz"))
        m_indexStress[2] = m_particles->getParamId("s_zz");
    else
        m_indexStress[2] = m_particles->registerParameter("s_zz");

    if(m_particles->hasParameter("s_xy"))
        m_indexStress[3] = m_particles->getParamId("s_xy");
    else
        m_indexStress[3] = m_particles->registerParameter("s_xy");

    if(m_particles->hasParameter("s_xz"))
        m_indexStress[4] = m_particles->getParamId("s_xz");
    else
        m_indexStress[4] = m_particles->registerParameter("s_xz");

    if(m_particles->hasParameter("s_yz"))
        m_indexStress[5] = m_particles->getParamId("s_yz");
    else
        m_indexStress[5] = m_particles->registerParameter("s_yz");

    m_state = false;
    m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
    m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepTwo(const pair<int, int> &id_col)
{
    // First calculating the total stress on an material point. The stress is averaged
    // all its bonds, and the stress state on a bond is the mean of the
    // stress at each material point.

    const int id_i = id_col.first;
    const int col_i = id_col.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    /*
    if(m_dim == 2)
    {
        arma::mat S_i(2, 2);
        arma::mat S(2, 2);
        S_i(0, 0) = (*m_data)(col_i, m_indexStress[0]);
        S_i(1, 1) = (*m_data)(col_i, m_indexStress[1]);
        S_i(0, 1) = (*m_data)(col_i, m_indexStress[3]);
        S_i(1, 0) = S_i(0, 1);
        arma::vec eigval;

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_pIds)[id_j];

            if((*m_data)(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(S_i(0, 0) + (*m_data)(j, m_indexStress[0]));
            S(1, 1) = 0.5*(S_i(1, 1) + (*m_data)(j, m_indexStress[1]));
            S(0, 1) = 0.5*(S_i(0, 1) + (*m_data)(j, m_indexStress[3]));
//            S(0, 0) = max(S_i(0, 0), (*m_data)(j, m_indexStress[0]));
//            S(1, 1) = max(S_i(1, 1), (*m_data)(j, m_indexStress[1]));
//            S(0, 1) = max(S_i(0, 1), (*m_data)(j, m_indexStress[3]));
            S(1, 0) = S(0, 1);
            const double first = 0.5*(S(0, 0) + S(1, 1));
            const double second = sqrt(0.25*(S(0, 0) - S(1, 1))*(S(0, 0) - S(1, 1)) + S(0, 1)*S(0, 1));
            const double s1 = first + second;
            const double s2 = first - second;
            const double p_1 = max(s1, s2);
            const double p_2 = min(s1, s2);

//            arma::eig_sym(eigval, S);
//            const double p_1 = eigval(1);
//            const double p_2 = eigval(0);

//            if(m_d*p_1 - p_2 - m_C > 0)
//            {
//                con.second[m_indexConnected] = 0;
//                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
//            }
//            else if(p_1 > m_T)
            if(p_1 > m_T)
            {
//                con.second[m_indexConnected] = 0;
                if(p_1 > m_maxStress)
                {
                    m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
                    m_maxStress = p_1;
                }
            }
        }
    }
    */
    arma::vec eigval(m_dim);
    arma::mat S_i(m_dim, m_dim);
    arma::mat S(m_dim, m_dim);

    if(m_dim == 2)
    {
        S_i(0, 0) = (*m_data)(col_i, m_indexStress[0]);
        S_i(1, 1) = (*m_data)(col_i, m_indexStress[1]);
        S_i(0, 1) = (*m_data)(col_i, m_indexStress[3]);
        S_i(1, 0) = S_i(0, 1);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int j = (*m_pIds)[id_j];

            if((*m_data)(j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(S_i(0, 0) + (*m_data)(j, m_indexStress[0]));
            S(1, 1) = 0.5*(S_i(1, 1) + (*m_data)(j, m_indexStress[1]));
            S(0, 1) = 0.5*(S_i(0, 1) + (*m_data)(j, m_indexStress[3]));
            S(1, 0) = S(0, 1);
#if 0
            arma::eig_sym(eigval, S);
            const double p_1 = eigval(1);
            const double p_2 = eigval(0);
#else
            const double first = 0.5*(S(0, 0) + S(1, 1));
            const double second = sqrt(0.25*(S(0, 0) - S(1, 1))*(S(0, 0) - S(1, 1)) + S(0, 1)*S(0, 1));
            const double s1 = first + second;
            const double s2 = first - second;
            const double p_1 = max(s1, s2);
            const double p_2 = min(s1, s2);
#endif

            if(m_d*p_1 - p_2 - m_C > 0)
            {
                con.second[m_indexConnected] = 0;
                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
            }
            else if(p_1 > m_T)
            {
                con.second[m_indexConnected] = 0;
                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
            }
        }
    }
    else if(m_dim == 3)
    {
        S_i(0, 0) = (*m_data)(col_i, m_indexStress[0]);
        S_i(1, 1) = (*m_data)(col_i, m_indexStress[1]);
        S_i(2, 2) = (*m_data)(col_i, m_indexStress[2]);
        S_i(0, 1) = (*m_data)(col_i, m_indexStress[3]);
        S_i(1, 0) = S_i(0, 1);
        S_i(0, 2) = (*m_data)(col_i, m_indexStress[4]);
        S_i(2, 0) = S_i(0, 2);
        S_i(1, 2) = (*m_data)(col_i, m_indexStress[5]);
        S_i(2, 1) = S_i(1, 2);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = (*m_pIds)[id_j];

            if((*m_data)(col_j, m_indexUnbreakable) >= 1)
                continue;

            if(con.second[m_indexConnected] <= 0.5)
                continue;

            S(0, 0) = 0.5*(S_i(0, 0) + (*m_data)(col_j, m_indexStress[0]));
            S(1, 1) = 0.5*(S_i(1, 1) + (*m_data)(col_j, m_indexStress[1]));
            S(2, 2) = 0.5*(S_i(2, 2) + (*m_data)(col_j, m_indexStress[2]));
            S(0, 1) = 0.5*(S_i(0, 1) + (*m_data)(col_j, m_indexStress[3]));
            S(1, 0) = S(0, 1);
            S(0, 2) = 0.5*(S_i(0, 2) + (*m_data)(col_j, m_indexStress[4]));
            S(2, 0) = S(0, 2);
            S(1, 2) = 0.5*(S_i(1, 2) + (*m_data)(col_j, m_indexStress[5]));
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
                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
            }
            else if(p_1 > m_T)
            {
                con.second[m_indexConnected] = 0;
                m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
            }
        }
    }
    //--------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepOne(const pair<int, int> &pIdcol)
//void ADRmohrCoulombFracture::evaluateStepTwo(const pair<int, int> &pIdcol)
{
    for(int s=0; s<6; s++)
    {
        (*m_data)(pIdcol.second, m_indexStress[s]) = 0;
    }

    for(Force *force: m_forces)
    {
        force->calculateStress(pIdcol, m_indexStress);
    }
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepTwo()
{
    if(m_maxPId.first != -1)
    {
        m_state = true;

        auto &con = m_maxPId.second;
        con->second[m_indexConnected] = 0;
    }
    else
    {
        m_state = false;
    }

    m_maxPId = std::pair<int, std::pair<int, vector<double>> *>(-1, nullptr);
    m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::addForce(Force *force)
{
    m_forces.push_back(force);
}
//------------------------------------------------------------------------------
}

