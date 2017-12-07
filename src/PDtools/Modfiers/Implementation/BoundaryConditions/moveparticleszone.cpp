#include "moveparticleszone.h"

#include "PDtools/Particles/pd_particles.h"
namespace PDtools {
//------------------------------------------------------------------------------
MoveParticlesZone::MoveParticlesZone(double velAmplitude, vec velocityDirection,
                                     vector<double> initialArea, double dt,
                                     bool isStatic, double delta)
    : m_delta(delta), m_velAmplitude(velAmplitude), m_initialArea(initialArea),
      m_dt(dt), m_isStatic(isStatic), m_time(dt),
      m_usingUnbreakableBorder(true),
      m_velocityDirection(arma::normalise(velocityDirection)) {}
//------------------------------------------------------------------------------
void MoveParticlesZone::registerParticleParameters() {
  m_particles->registerParameter("unbreakable");
  m_ghostParameters = {"unbreakable"};
}
//------------------------------------------------------------------------------
void MoveParticlesZone::evaluateStepOne() {
  const unordered_map<int, int> &idToCol = m_particles->idToCol();
  //  mat &r = m_particles->r();
  const double v_amp = m_time * m_velAmplitude;
  const vec dr = v_amp * m_velocityDirection;
  arma::mat &v = m_particles->v();
  arma::mat &F = m_particles->F();
  arma::mat &Fold = m_particles->Fold();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol.at(id);
    for (int d = 0; d < m_dim; d++) {
      //            r(i, d) += dr(d);
      v(i, d) = (1. - m_vd_abs(d)) * v(i, d) +
                m_velocityDirection(d) * m_velAmplitude;
      F(i, d) = (1. - m_vd_abs(d)) * F(i, d);
      Fold(i, d) = (1. - m_vd_abs(d)) * Fold(i, d);
    }
  }
}
//------------------------------------------------------------------------------
void MoveParticlesZone::staticEvaluation() {
  const unordered_map<int, int> &idToCol = m_particles->idToCol();
  arma::mat &v = m_particles->v();
  arma::mat &F = m_particles->F();
  arma::mat &Fold = m_particles->Fold();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol.at(id);
    for (int d = 0; d < m_dim; d++) {
      v(i, d) = (1 - m_vd_abs(d)) * v(i, d);
      F(i, d) = (1 - m_vd_abs(d)) * F(i, d);
      Fold(i, d) = (1 - m_vd_abs(d)) * Fold(i, d);
    }
  }
}
//------------------------------------------------------------------------------
void MoveParticlesZone::initialize() {
  for (int d = 0; d < m_dim; d++)
    m_vd_abs(d) = fabs(m_velocityDirection(d));

  // Selecting particles
  const ivec &colToId = m_particles->colToId();
  const arma::mat &r = m_particles->r();
  arma::mat &data = m_particles->data();
  arma::imat &isStatic = m_particles->isStatic();
  const int unbreakablePos = m_particles->registerParameter("unbreakable");

  //    double uRadius = 6*m_delta;
  double uRadius = 4.5 * m_delta;
  int n = 0;

  for (unsigned int i = 0; i < m_particles->nParticles(); i++) {
    bool inArea = true;
    bool unbreakable = false;
    if (m_usingUnbreakableBorder)
      unbreakable = true;

    for (int d = 0; d < m_dim; d++) {
      if (m_initialArea[2 * d] <= r(i, d) &&
          r(i, d) < m_initialArea[2 * d + 1]) {
        // in area
      } else
        inArea = false;

      if (m_usingUnbreakableBorder) {
        if (m_initialArea[2 * d] - uRadius <= r(i, d) &&
            r(i, d) <=
                m_initialArea[2 * d + 1] + uRadius) { // in unbreakable area
        } else {
          unbreakable = false;
        }
      }
    }

    if (inArea) {
      if (m_isStatic)
        isStatic(i) = 1;

      const int id = colToId(i);
      m_localParticleIds.push_back(id);
      data(i, unbreakablePos) = 1;
    }
    if (unbreakable) {
      data(i, unbreakablePos) = 1;
      n++;
    }
  }

  //    cout << "done: " << m_velAmplitude  << endl;
}
//------------------------------------------------------------------------------
}
