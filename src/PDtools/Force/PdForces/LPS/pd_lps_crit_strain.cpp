#include "pd_lps_crit_strain.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PD_LPS_CRIT_STRAIN::PD_LPS_CRIT_STRAIN(PD_Particles &particles, double c, double stretchCrit, double shearCrit, bool planeStress, bool analyticalM):
    PD_LPS(particles, planeStress, analyticalM),
    m_stretchCrit(stretchCrit),
    m_shearCrit(shearCrit),
    m_dampCoeff(c)
{
    m_hasStepTwoModifier = true;
    m_iUnbreakable = particles.registerParameter("unbreakable");
    m_particles.registerParameter("damage");
}
//------------------------------------------------------------------------------
void PD_LPS_CRIT_STRAIN::evaluateStepTwo(int id_i, int i)
{
    if(m_data(i, m_iUnbreakable) >= 1)
        return;

    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);
    bool broken = false;
//    const double theta_i = m_data(i, m_iTheta);

    for(auto &con:PDconnections) {
        const int id_j = con.first;
        const int j = m_idToCol[id_j];

        if(m_data(j, m_iUnbreakable) >= 1)
            continue;

        if(con.second[m_iConnected] <= 0.5)
            continue;

        const double s = con.second[m_iStretch];
        if(s > m_stretchCrit) {
            m_data(i, m_indexBrokenNow) = 1;
            con.second[m_iConnected] = 0;
            m_continueState = true;
            broken = true;
        }
        //        const double theta_j = m_data(j, m_iTheta);
        //        const double theta = 0.5*(theta_i + theta_j);
        //        const double s_i = theta/3.;
        //        const double s_d = s - s_i;
        //        const double s_d = max(s - theta_i/3., s - theta_j/3.);
        //        else if(s_d > m_shearCrit) {
        //            m_data(i, m_indexBrokenNow) = 1;
        //            con.second[m_iConnected] = 0;
        //            m_continueState = true;
        //            broken = true;
        //        }
    }

    if(broken) {
        updateWeightedVolume(id_i, i);
//        m_data(i, m_indexBrokenNow) = 0;
    }
}
//------------------------------------------------------------------------------
}
