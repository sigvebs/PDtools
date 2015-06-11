#include "computestress.h"
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
ComputeStress::ComputeStress(PD_Particles &particles, vector<Force *> &forces):
    ComputeProperty(particles),
    m_forces(forces),
    m_data(m_particles.data())
{
    if(m_particles.hasParameter("s_xx"))
        m_indexStress[0] = m_particles.getParamId("s_xx");
    else
        m_indexStress[0] = m_particles.registerParameter("s_xx");

    if(m_particles.hasParameter("s_yy"))
        m_indexStress[1] = m_particles.getParamId("s_yy");
    else
        m_indexStress[1] = m_particles.registerParameter("s_yy");

    if(m_particles.hasParameter("s_zz"))
        m_indexStress[2] = m_particles.getParamId("s_zz");
    else
        m_indexStress[2] = m_particles.registerParameter("s_zz");

    if(m_particles.hasParameter("s_xy"))
        m_indexStress[3] = m_particles.getParamId("s_xy");
    else
        m_indexStress[3] = m_particles.registerParameter("s_xy");

    if(m_particles.hasParameter("s_xz"))
        m_indexStress[4] = m_particles.getParamId("s_xz");
    else
        m_indexStress[4] = m_particles.registerParameter("s_xz");

    if(m_particles.hasParameter("s_yz"))
        m_indexStress[5] = m_particles.getParamId("s_yz");
    else
        m_indexStress[5] = m_particles.registerParameter("s_yz");

//    m_indexStress[1] = m_particles.registerParameter("s_yy");
//    m_indexStress[2] = m_particles.registerParameter("s_zz");
//    m_indexStress[3] = m_particles.registerParameter("s_xy");
//    m_indexStress[4] = m_particles.registerParameter("s_xz");
//    m_indexStress[5] = m_particles.registerParameter("s_yz");
}
//------------------------------------------------------------------------------
ComputeStress::~ComputeStress()
{

}
//------------------------------------------------------------------------------
void ComputeStress::update(const pair<int, int> &pIdcol)
{
    for(int s=0; s<6; s++)
    {
        m_data(pIdcol.second, m_indexStress[s]) = 0;
    }

    for(Force *force: m_forces)
    {
        force->calculateStress(pIdcol, m_indexStress);
    }
}
//------------------------------------------------------------------------------
}


