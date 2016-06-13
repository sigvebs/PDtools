#include "moveparticletype.h"

#include "PDtools/Particles/pd_particles.h"
namespace PDtools
{
//------------------------------------------------------------------------------
MoveParticleGroup::MoveParticleGroup(int groupId, double velAmplitude, vec velocityDirection,
                                   double dt, bool isStatic):
    m_groupId(groupId),
    m_velAmplitude(velAmplitude),
    m_dt(dt),
    m_isStatic(isStatic),
    m_time(dt),
    m_velocityDirection(arma::normalise(velocityDirection))
{
}
//------------------------------------------------------------------------------
void MoveParticleGroup::registerParticleParameters()
{
    m_particles->registerParameter("unbreakable");
    m_ghostParameters = {"unbreakable"};
}
//------------------------------------------------------------------------------
void MoveParticleGroup::evaluateStepOne()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    mat & r = m_particles->r();
    const vec dr = m_time*m_velAmplitude*m_velocityDirection;
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i =  idToCol.at(id);
        for(int d=0; d<m_dim; d++)
        {
            r(i, d) += dr(d);
            v(i, d) = (1. - m_vd(d))*v(i, d);
            F(i, d) = (1. - m_vd(d))*F(i, d);
            Fold(i, d) = (1. - m_vd(d))*Fold(i, d);
        }
    }
}
//------------------------------------------------------------------------------
void MoveParticleGroup::staticEvaluation()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & v = m_particles->v();
    arma::mat & F = m_particles->F();
    arma::mat & Fold = m_particles->Fold();

    for(const int &id:m_localParticleIds)
    {
        const int i = idToCol.at(id);
        for(int d=0; d<m_dim; d++)
        {
            v(i, d) = (1 - m_vd(d))*v(i, d);
            F(i, d) = (1 - m_vd(d))*F(i, d);
            Fold(i, d) = (1 - m_vd(d))*Fold(i, d);
        }
    }
}
//------------------------------------------------------------------------------
void MoveParticleGroup::initialize()
{
    for(int d=0; d<m_dim; d++)
        m_vd(d) = fabs(m_velocityDirection(d));

    // Selecting particles
    const ivec &colToId = m_particles->colToId();
    arma::mat & data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();
    const int iType = m_particles->getParamId("groupId");
    const int iUnbreakable = m_particles->registerParameter("unbreakable");

    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        const int type = data(i, iType);

        if(type == m_groupId)
        {
            if(m_isStatic)
                isStatic(i) = 1;

            const int id = colToId(i);
            m_localParticleIds.push_back(id);
            data(i, iUnbreakable) = 1;
        }
    }
}
//------------------------------------------------------------------------------
}
