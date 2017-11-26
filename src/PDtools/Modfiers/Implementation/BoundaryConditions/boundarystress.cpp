#include "boundarystress.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
BoundaryStress::BoundaryStress(double appliedStress, double stressOrientation,
                               pair<double, double> boundary, int boundaryOrientation,
                               int steps, double delta)
{
    m_appliedStress = appliedStress;
    m_stressOritentation = stressOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_inceremental = false;
    m_incrementalStress = 0;
    m_delta = delta;
    if(!m_inceremental) {
        m_incrementalStress = m_appliedStress;
    }
}
//------------------------------------------------------------------------------
void BoundaryStress::evaluateStepOne()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & F = m_particles->F();
    arma::mat & data = m_particles->data();

    for(const int &id:m_localParticleIds) {
        const int i = idToCol.at(id);
        const double radius = data(i, m_indexRadius);
        const double volume = data(i, m_indexVolume);
        const double area = M_PI*radius*radius;
        const double f = m_incrementalStress*area/volume;
        F(i, m_stressOritentation) += f;
    }
}
//------------------------------------------------------------------------------
void BoundaryStress::evaluateStepTwo()
{
    arma::mat & b = m_particles->b();
    b.zeros();
    if(m_inceremental) {
        m_incrementalStress += m_appliedStress;
    }
//    evaluateStepOne();
}
//------------------------------------------------------------------------------
void BoundaryStress::staticEvaluation()
{
    const unordered_map<int, int> &idToCol = m_particles->idToCol();
    arma::mat & F = m_particles->F();
    arma::mat & data = m_particles->data();

    for(const int &id:m_localParticleIds) {
        const int i = idToCol.at(id);

        const double radius = data(i, m_indexRadius);
        const double area = M_PI*radius*radius;
        const double volume = data(i, m_indexVolume);
        const double f = m_incrementalStress*area/volume;
        F(i, m_stressOritentation) += f;
    }
}
//------------------------------------------------------------------------------
void BoundaryStress::initialize()
{
    const ivec &colToId = m_particles->colToId();
    m_indexRadius = m_particles->registerParameter("radius");
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();

    m_indexVolume = m_particles->getParamId("volume");
    double volume = 0;
    int unbreakablePos;

    if(m_particles->hasParameter("unbreakable")) {
        unbreakablePos = m_particles->getParamId("unbreakable");
    } else {
        unbreakablePos = m_particles->registerParameter("unbreakable");
    }

    double uRadius = 0.15*(m_boundary.second - m_boundary.first);

    // Selecting particles
    for(unsigned int i=0; i<m_particles->nParticles(); i++) {
        const double pos = r(i, m_boundaryOrientation);

        if(m_boundary.first <= pos && pos < m_boundary.second) {
            const int id = colToId(i);
            m_localParticleIds.push_back(id);
            volume += data(i, m_indexVolume);
            data(i, unbreakablePos) = 1;
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second) {
            data(i, unbreakablePos) = 1;
        }
    }
}
//------------------------------------------------------------------------------
}
