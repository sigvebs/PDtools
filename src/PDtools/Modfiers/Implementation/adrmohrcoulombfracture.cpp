#include "adrmohrcoulombfracture.h"

#include "PDtools/Force/force.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ADRmohrCoulombFracture::ADRmohrCoulombFracture(double mu, double C, double T):
    m_C(C), m_T(T)
{
    m_d = pow(sqrt(1 + mu*mu) + mu, 2);
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

    m_dim = 2;
    m_state = false;
    m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
    m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepOne(const pair<int, int> &pIdcol)
{
    // First calculating the total stress on an material point. The stress is averaged
    // all its bonds, and the stress state on a bond is the mean of the
    // stress at each material point.

    int id_i = pIdcol.first;
    int col_i = pIdcol.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
    vector<pair<int, vector<double>> *> removeParticles;


    if(m_dim == 2)
    {
        double s_xx_i = (*m_data)(col_i, m_indexStress[0]);
        double s_yy_i = (*m_data)(col_i, m_indexStress[1]);
        double s_xy_i = (*m_data)(col_i, m_indexStress[3]);

        for(auto &con:PDconnections)
        {
            int id_j = con.first;
            int col_j = (*m_pIds)[id_j];

            if((*m_data)(col_j, m_indexUnbreakable) >= 1)
                continue;

            double s_xx_j = (*m_data)(col_j, m_indexStress[0]);
            double s_yy_j = (*m_data)(col_j, m_indexStress[1]);
            double s_xy_j = (*m_data)(col_j, m_indexStress[3]);

            double sx = 0.5*(s_xx_i + s_xx_j);
            double sy = 0.5*(s_yy_i + s_yy_j);
            double s_xy = 0.5*(s_xy_i + s_xy_j);

            double first = 0.5*(sx + sy);
            double second = sqrt(0.25*(sx - sy)*(sx - sy) + s_xy*s_xy);
            double s1 = first + second;
            double s2 = first - second;

            double pStressOne = 0;
            double pStressTwo = 0;

            if(s1 > s2)
            {
                pStressOne = s1;
                pStressTwo = s2;
            }
            else
            {
                pStressOne = s2;
                pStressTwo = s1;
            }

            if(pStressOne > m_T)
            {
                if(pStressOne > m_maxStress)
                {
                    m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
                    m_maxStress = pStressOne;
                }
            }
//            if(m_d*pStressOne - pStressTwo  > -m_C)
//            {
//                removeParticles.push_back(&con);
//            }
//            else if(pStressOne > m_T)
//            {
//                removeParticles.push_back(&con);
//            }
        }
    }
    //--------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::evaluateStepTwo(const pair<int, int> &pIdcol)
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
        int pId = m_maxPId.first;
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
        PDconnections.erase(remove(begin(PDconnections), end(PDconnections), *m_maxPId.second),
                             end(PDconnections) );
    }
    else
    {
        m_state = false;
    }

    m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
    m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
void ADRmohrCoulombFracture::addForce(Force *force)
{
    m_forces.push_back(force);
}
//------------------------------------------------------------------------------
}

