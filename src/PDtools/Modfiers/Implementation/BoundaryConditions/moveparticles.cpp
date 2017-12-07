#include "moveparticles.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
MoveParticles::MoveParticles(double velAmplitude, double velOrientation,
                             pair<double, double> boundary,
                             int boundaryOrientation, double dt,
                             bool isStatic) {
  m_velAmplitude = velAmplitude;
  m_velOritentation = velOrientation;
  m_boundary = boundary;
  m_boundaryOrientation = boundaryOrientation;
  m_dt = dt;
  m_isStatic = isStatic;
  m_time = dt;
  m_usingUnbreakableBorder = true;
}
//------------------------------------------------------------------------------
void MoveParticles::registerParticleParameters() {
  m_particles->registerParameter("unbreakable");
  m_ghostParameters = {"unbreakable"};

  m_particles->setNeedGhostVelocity(1);
  m_particles->setNeedGhostR0(1);
}
//------------------------------------------------------------------------------
void MoveParticles::evaluateStepOne() {
  const unordered_map<int, int> &idToCol = m_particles->idToCol();
  mat &r = m_particles->r();
  const double dr = m_time * m_velAmplitude;
  arma::mat &v = m_particles->v();
  arma::mat &F = m_particles->F();
  arma::mat &Fold = m_particles->Fold();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol.at(id);
    r(i, m_velOritentation) += dr;
    v(i, m_velOritentation) = 0.0;
    F(i, m_velOritentation) = 0.0;
    Fold(i, m_velOritentation) = 0.0;
  }
}
//------------------------------------------------------------------------------
void MoveParticles::initialize() {
  // Selecting particles
  const ivec &colToId = m_particles->colToId();
  const arma::mat &r = m_particles->r();
  arma::mat &data = m_particles->data();
  arma::imat &isStatic = m_particles->isStatic();
  const int unbreakablePos = m_particles->registerParameter("unbreakable");

  double uRadius = 0.15 * (m_boundary.second - m_boundary.first);

  // Selecting particles
  for (unsigned int i = 0; i < m_particles->nParticles(); i++) {
    const double pos = r(i, m_boundaryOrientation);

    if (m_boundary.first <= pos && pos < m_boundary.second) {
      const int id = colToId(i);

      if (m_isStatic)
        isStatic(i) = 1;

      m_localParticleIds.push_back(id);
      data(i, unbreakablePos) = 1;
    }
    if (m_usingUnbreakableBorder) {
      if (m_boundary.first - uRadius <= pos &&
          pos < uRadius + m_boundary.second) {
        data(i, unbreakablePos) = 1;
      }
    }
  }
}
//------------------------------------------------------------------------------
void MoveParticles::staticEvaluation() {
  const unordered_map<int, int> &idToCol = m_particles->idToCol();
  arma::mat &v = m_particles->v();
  arma::mat &F = m_particles->F();
  arma::mat &Fold = m_particles->Fold();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol.at(id);

    v(i, m_velOritentation) = 0.0;
    F(i, m_velOritentation) = 0.0;
    Fold(i, m_velOritentation) = 0.0;
  }
}
//------------------------------------------------------------------------------
}
