#include "bondenergyfracture.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
BondEnergyFracture::BondEnergyFracture(double delta, double G0, int dim,
                                       vector<Force*> &forces, double h):
    m_forces(forces)
{
    m_G = G0;
    if(dim == 3)
    {
        m_wc = 4.*G0 / (M_PI*pow(delta, 4));
    }else
    {
        m_wc = 9.*G0 / (4.*h*pow(delta, 3));
    }
}
//------------------------------------------------------------------------------
BondEnergyFracture::~BondEnergyFracture()
{

}
//------------------------------------------------------------------------------
void BondEnergyFracture::initialize()
{
    m_indexUnbreakable =  m_particles->getParamId("unbreakable");
    m_indexStretch = m_particles->getPdParamId("stretch");
    m_indexDr0 = m_particles->getPdParamId("dr0");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_indexForceScaling = m_particles->getPdParamId("forceScalingBond");
    m_indexMicromodulus = m_particles->getParamId("micromodulus");
    m_pIds = &m_particles->pIds();
    m_data = &m_particles->data();

    m_state = false;
    m_broken = false;
    for(Force *force:m_forces)
    {
        force->updateState();
    }
}
//------------------------------------------------------------------------------
void BondEnergyFracture::evaluateStepOne(const pair<int, int> &pIdcol)
{
    int id_i = pIdcol.first;
    int col_i = pIdcol.second;

    if((*m_data)(col_i, m_indexUnbreakable) >= 1)
        return;

    const double c_i = (*m_data)(col_i, m_indexMicromodulus);

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        int id_j = con.first;
        int col_j = (*m_pIds)[id_j];

        if((*m_data)(col_j, m_indexUnbreakable) >= 1)
            continue;
        if(con.second[m_indexConnected] <= 0.5)
            continue;

        const double c_j = (*m_data)(col_j, m_indexMicromodulus);
        const double g_ij = con.second[m_indexForceScaling];
        const double c = 0.5*(c_i + c_j)*g_ij;
        const double s = con.second[m_indexStretch];
        const double dr0 = con.second[m_indexDr0];
//        const double w = 0.5*c*s*s*dr0;
        double w = 0;
//        for(Force *force:m_forces)
//        {
//            force->updateState();
////            w += force->calculateBondEnergy(pIdcol, con);
//        }
        if(w > m_wc)
        {
            con.second[m_indexConnected] = 0;
            m_broken = true;
        }

//        double sc = sqrt(m_G/dr0);
//        if(s > sc)
//        {
//            con.second[m_indexConnected] = 0;
//            m_broken = true;
//        }
    }
}
//------------------------------------------------------------------------------
void BondEnergyFracture::evaluateStepTwo()
{
    if(m_broken)
    {
        m_state = true;
    }
    else
    {
        m_state = false;
    }
    m_broken = false;
}
//------------------------------------------------------------------------------
}

