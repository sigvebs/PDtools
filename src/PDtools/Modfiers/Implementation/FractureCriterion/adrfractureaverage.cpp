#include "adrfractureaverage.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
ADRfractureAverage::ADRfractureAverage(double alpha) : m_alpha(alpha) {
  m_hasStepOne = true;
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void ADRfractureAverage::initialize() {
  m_indexS0 = m_particles->getParamId("s0");
  m_indexUnbreakable = m_particles->getParamId("unbreakable");
  m_indexStretch = m_particles->getPdParamId("stretch");
  m_indexS00 = m_particles->registerPdParameter("s00");
  m_indexS_avg = m_particles->registerParameter("s_avg");
  m_idToCol = &m_particles->getIdToCol_v();
  m_data = &m_particles->data();
  const ivec &colToId = m_particles->colToId();

  // Setting the initial max stretch between two particles
  for (unsigned int i = 0; i < m_particles->nParticles(); i++) {
    const int pId = colToId(i);
    const double s0_i = (*m_data)(i, m_indexS0);

    vector<pair<int, vector<double>>> &PDconnections =
        m_particles->pdConnections(pId);
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int col_j = (*m_idToCol)[id_j];
      const double s0_j = (*m_data)(col_j, m_indexS0);
      con.second[m_indexS00] = 0.5 * (s0_i + s0_j);
    }
  }

  m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
  m_maxStretch = std::numeric_limits<double>::min();
  m_state = false;
}
//------------------------------------------------------------------------------
void ADRfractureAverage::evaluateStepOne(const int id_i, const int i) {
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;

  const double s0_i = (*m_data)(i, m_indexS0);
  const double s_i = (*m_data)(i, m_indexS_avg);
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = (*m_idToCol)[id_j];
    if ((*m_data)(j, m_indexUnbreakable) >= 1)
      continue;

    const double s0_j = (*m_data)(j, m_indexS0);
    const double s_j = (*m_data)(j, m_indexS_avg);
    const double s = 0.5 * (s_i + s_j);
    const double s0 = std::min(s0_i, s0_j);

    if (s > s0) {
      if (s > m_maxStretch) {
        m_maxPId = pair<int, pair<int, vector<double>> *>(id_i, &con);
        m_maxStretch = s;
      }
    }
  }
}
//------------------------------------------------------------------------------
void ADRfractureAverage::evaluateStepTwo(const pair<int, int> &pIdcol) {
  const int pId = pIdcol.first;
  const int col_i = pIdcol.second;

  if ((*m_data)(col_i, m_indexUnbreakable) >= 1)
    return;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(pId);

  double s_avg = 0;
  for (auto &con : PDconnections) {
    s_avg += con.second[m_indexStretch];
  }

  (*m_data)(col_i, m_indexS_avg) = s_avg / PDconnections.size();
}
//------------------------------------------------------------------------------
void ADRfractureAverage::evaluateStepTwo() {
  if (m_maxPId.first != -1) {
    m_state = true;
    int pId = m_maxPId.first;
    vector<pair<int, vector<double>>> &PDconnections =
        m_particles->pdConnections(pId);
    PDconnections.erase(
        remove(begin(PDconnections), end(PDconnections), *m_maxPId.second),
        end(PDconnections));
  } else {
    m_state = false;
  }

  m_maxPId = pair<int, pair<int, vector<double>> *>(-1, nullptr);
  m_maxStretch = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
}
