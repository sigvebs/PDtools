#include "boundaryforce.h"
#include "PDtools/Particles/pd_particles.h"

#if USE_MPI
#include <mpi.h>
#endif

namespace PDtools {
//------------------------------------------------------------------------------
boundaryForce::boundaryForce(double appliedForce, double forceOrientation,
                             pair<double, double> boundary,
                             int boundaryOrientation, int steps, double delta,
                             int incremental, int scale) {
  (void) steps;
  m_forceDensity = appliedForce;
  m_forceOritentation = forceOrientation;
  m_boundary = boundary;
  m_boundaryOrientation = boundaryOrientation;
  m_inceremental = incremental;
  m_incrementalForce = m_forceDensity;
  m_delta = delta;
  m_scale = scale;
  if (!m_inceremental) {
    m_incrementalForce = m_forceDensity;
  }
}
//------------------------------------------------------------------------------
void boundaryForce::registerParticleParameters() {
  m_particles->registerParameter("unbreakable");
  m_ghostParameters = {"unbreakable"};

  m_particles->setNeedGhostVelocity(1);
  m_particles->setNeedGhostR0(1);
}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepOne() {
  const ivec &idToCol = m_particles->getIdToCol_v();
  arma::mat &F = m_particles->F();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol[id];
    const double f = m_incrementalForce;
    F(i, m_forceOritentation) += f;
  }
}
//------------------------------------------------------------------------------
void boundaryForce::initialize() {
  const ivec &colToId = m_particles->colToId();
  m_indexRadius = m_particles->registerParameter("radius");
  const arma::mat &r = m_particles->r();
  arma::mat &data = m_particles->data();

  m_indexVolume = m_particles->getParamId("volume");
  double volume = 0;
  int unbreakablePos;

  if (m_particles->hasParameter("unbreakable")) {
    unbreakablePos = m_particles->getParamId("unbreakable");
  } else {
    unbreakablePos = m_particles->registerParameter("unbreakable");
  }

  double uRadius = 0.15 * (m_boundary.second - m_boundary.first);

  // Selecting particles
  for (unsigned int i = 0; i < m_particles->nParticles(); i++) {
    const double pos = r(i, m_boundaryOrientation);

    if (m_boundary.first <= pos && pos < m_boundary.second) {
      const int id = colToId(i);
      m_localParticleIds.push_back(id);
      volume += data(i, m_indexVolume);
      //            data(i, unbreakablePos) = 1;
    }

    if (m_boundary.first - uRadius <= pos &&
        pos < uRadius + m_boundary.second) {
      //            data(i, unbreakablePos) = 1;
    }
  }
  if (m_scale) {
#if USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    m_forceDensity /= volume;
    m_incrementalForce = m_forceDensity;
    if (!m_inceremental) {
      m_incrementalForce = m_forceDensity;
    }
  }
}
//------------------------------------------------------------------------------
void boundaryForce::evaluateStepTwo() {
  if (m_inceremental) {
    m_incrementalForce += m_forceDensity;
  }
  //    evaluateStepOne();
}
//------------------------------------------------------------------------------
void boundaryForce::staticEvaluation() {
  const ivec &idToCol = m_particles->getIdToCol_v();
  arma::mat &F = m_particles->F();

  for (const int &id : m_localParticleIds) {
    const int i = idToCol[id];
    const double f = m_incrementalForce;
    F(i, m_forceOritentation) += f;
  }
  //    boundaryForce::evaluateStepTwo();
}
//------------------------------------------------------------------------------
}
