#include "mohrcoulombmaxfractureweighted.h"
#include "PDtools/Grid/grid.h"
#include "PDtools/Particles/pd_particles.h"

#include <stdlib.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace PDtools {
//------------------------------------------------------------------------------
MohrCoulombMaxFractureWeighted::MohrCoulombMaxFractureWeighted(
    double mu, double C, double T, double Wc, double Bc)
    : m_C(C), m_T(T) {
  m_phi = mu * M_PI / 180.;
  m_d = tan(m_phi);
  m_neededProperties = {pair<string, int>("stress", 1)};

  m_Wc = Wc;
  m_Bc = Bc;
  m_Wn = 1. - m_Wc;
  m_Bn = 1. - m_Bc;

  m_cosTheta = cos(M_PI / 2. + m_phi);
  m_sinTheta = sin(M_PI / 2. + m_phi);
  m_radiusScale = 1.15;

  m_hasStepOne = true;
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFractureWeighted::initialize() {
  srand(time(NULL));
  m_broken = false;
  m_brokenParticles.clear();
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFractureWeighted::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexRadius = m_particles->registerParameter("radius");
  m_indexCompute = m_particles->registerPdParameter("compute");
  m_indexBroken = m_particles->registerParameter("broken", 0);
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_particles->registerParameter("damage");
  m_idToCol = &m_particles->getIdToCol_v();

  switch (m_dim) {
  case 1:
    m_ghostParameters = {"s_xx"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    break;
  case 2:
    m_ghostParameters = {"s_xx", "s_yy", "s_xy"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    break;
  case 3:
    m_ghostParameters = {"s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStress[3] = m_particles->registerParameter("s_zz");
    m_indexStress[4] = m_particles->registerParameter("s_xz");
    m_indexStress[5] = m_particles->registerParameter("s_yz");
    break;
  }
  m_ghostParameters.push_back("broken");
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFractureWeighted::evaluateStepOne(const int id_i,
                                                     const int i) {
  mat &data = *m_data;

  if (data(i, m_indexUnbreakable) >= 1)
    return;

  const int broken_i = data(i, m_indexBroken);
  vector<int> broken_nodes;

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = (*m_idToCol)[id_j];

    if (data(j, m_indexUnbreakable) >= 1)
      continue;

    const int broken_j = data(j, m_indexBroken);

    if (con.second[m_indexConnected] <= 0.5)
      continue;

    if (broken_i > 0) {
      double sx, sy, sxy;
      sx = m_Wn * data(j, m_indexStress[0]) + m_Wc * data(i, m_indexStress[0]);
      sy = m_Wn * data(j, m_indexStress[1]) + m_Wc * data(i, m_indexStress[1]);
      sxy = m_Wn * data(j, m_indexStress[2]) + m_Wc * data(i, m_indexStress[2]);

      const double first = 0.5 * (sx + sy);
      const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

      const double p_2 = first + second;
      const double p_1 = first - second;

      const double shear = 0.5 * (p_1 - p_2) * m_sinTheta;
      const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * m_cosTheta;

      const double criticalShear = fabs(shear) - fabs(m_C - m_d * normal);
      const double criticalTensile = p_2 - m_T;

      if (criticalShear >= 0 && normal < 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
        m_broken = true;
      } else if (criticalTensile >= 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
        m_broken = true;
      }

      continue;
    }

    if (broken_j > 0) {
      broken_nodes.push_back(id_j);
      double sx, sy, sxy;
      sx = m_Wn * data(i, m_indexStress[0]) + m_Wc * data(j, m_indexStress[0]);
      sy = m_Wn * data(i, m_indexStress[1]) + m_Wc * data(j, m_indexStress[1]);
      sxy = m_Wn * data(i, m_indexStress[2]) + m_Wc * data(j, m_indexStress[2]);

      const double first = 0.5 * (sx + sy);
      const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

      const double p_2 = first + second; // max
      const double p_1 = first - second; // min

      const double shear = 0.5 * (p_1 - p_2) * m_sinTheta;
      const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * m_cosTheta;

      const double criticalShear = fabs(shear) - fabs(m_C - m_d * normal);
      const double criticalTensile = p_2 - m_T;

      if (criticalShear >= 0 && normal < 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
        m_broken = true;
      } else if (criticalTensile >= 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
        m_broken = true;
      }

      continue;
    }
  }

  if (!broken_nodes.empty()) {
    const mat &r = m_particles->r();
    double ip[m_dim];
    double ij[m_dim];
    double p[m_dim];

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol)[id_j];

      if (data(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      ij[0] = r(j, 0) - r(i, 0);
      ij[1] = r(j, 1) - r(i, 1);
      double ij2 = ij[0] * ij[0] + ij[1] * ij[1];

      for (int id_g : broken_nodes) {
        const int g = (*m_idToCol)[id_g];
        ip[0] = r(g, 0) - r(i, 0);
        ip[1] = r(g, 1) - r(i, 1);
        const double t = (ij[0] * ip[0] + ij[1] * ip[1]) / ij2;

        p[0] = r(i, 0) + ij[0] * t;
        p[1] = r(i, 1) + ij[1] * t;

        double dr2 = (p[0] - r(g, 0)) * (p[0] - r(g, 0)) +
                     (p[1] - r(g, 1)) * (p[1] - r(g, 1));
        const double radius_j = m_radiusScale * data(j, m_indexRadius);
        const double r2 = radius_j * radius_j;

        if (dr2 < r2) {
          if (j == g)
            continue;

          double ndr_i = 0;
          double ndr_j = 0;
          for (int d = 0; d < 2; d++) {
            //                        ndr_i += ij[d]*ip[d];
            //                        ndr_j -= ij[d]*(r(g, d) - r(j, d));
            ndr_i += (r(j, d) - r(i, d)) * (r(g, d) - r(i, d));
            ndr_j += (r(i, d) - r(j, d)) * (r(g, d) - r(j, d));
          }
          if (ndr_i < 0 || ndr_j < 0)
            continue;

          double sx, sy, sxy;
          sx = 0.5 * m_Bn *
                   (data(i, m_indexStress[0]) + data(j, m_indexStress[0])) +
               m_Bc * data(g, m_indexStress[0]);
          sy = 0.5 * m_Bn *
                   (data(i, m_indexStress[1]) + data(j, m_indexStress[1])) +
               m_Bc * data(g, m_indexStress[1]);
          sxy = 0.5 * m_Bn *
                    (data(i, m_indexStress[2]) + data(j, m_indexStress[2])) +
                m_Bc * data(g, m_indexStress[2]);

          const double first = 0.5 * (sx + sy);
          const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

          const double p_2 = first + second; // max
          const double p_1 = first - second; // min

          const double shear = 0.5 * (p_1 - p_2) * m_sinTheta;
          const double normal =
              0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * m_cosTheta;

          const double criticalShear = fabs(shear) - fabs(m_C - m_d * normal);
          const double criticalTensile = p_2 - m_T;
          int br = 0;
          if (criticalShear >= 0 && normal < 0) {
            con.second[m_indexConnected] = 0;
            data(i, m_indexBrokenNow) = 1;
            m_broken = true;
            br = 1;
          } else if (criticalTensile >= 0) {
            con.second[m_indexConnected] = 0;
            data(i, m_indexBrokenNow) = 1;
            m_broken = true;
            br = 1;
          }
        }
      }
    }
  }
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFractureWeighted::evaluateStepOnePost() {
#if USE_MPI
  exchangeBrokenParticlesMPI();
#endif
  const ivec &idToCol = m_particles->getIdToCol_v();
  const int nParticles = m_particles->nParticles();
  for (auto broken_idList : m_brokenParticles) {
    const int id_i = broken_idList.first;
    const vector<int> &broken = broken_idList.second;

#if USE_MPI
    // Check if the particle is a ghost particle
    if (idToCol[id_i] >= nParticles) {
      continue;
    }
#endif
    vector<pair<int, vector<double>>> &PDconnections =
        m_particles->pdConnections(id_i);

    for (auto &con : PDconnections) {
      for (int id : broken) {
        if (id == id_i) {
          con.second[m_indexConnected] = 0;
        }
      }
    }
  }
  m_brokenParticles.clear();
}
//------------------------------------------------------------------------------
// void MohrCoulombMaxFractureWeighted::evaluateStepOne(const int id_i, const
// int i)
void MohrCoulombMaxFractureWeighted::evaluateStepTwo(const int id_i,
                                                     const int i) {
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;

  mat &data = *m_data;
  mat &r = m_particles->r();

  if (m_dim == 2) {
    double sx, sy, sxy;
    sx = data(i, m_indexStress[0]);
    sy = data(i, m_indexStress[1]);
    sxy = data(i, m_indexStress[2]);

    const double first = 0.5 * (sx + sy);
    const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

    const double p_2 = first + second; // max
    const double p_1 = first - second; // min

    const double shear = 0.5 * (p_1 - p_2) * m_sinTheta;
    const double normal = 0.5 * (p_1 + p_2) + 0.5 * (p_1 - p_2) * m_cosTheta;

    const double criticalShear = fabs(shear) - fabs(m_C - m_d * normal);
    const double criticalTensile = p_2 - m_T;

    //        if(criticalShear >= 0  && normal < 0)
    if (criticalShear >= 0) {
      data(i, m_indexBroken) = 2;
      m_broken = true;
    } else if (criticalTensile >= 0) {
      data(i, m_indexBroken) = 1;
      m_broken = true;
    } else {
      data(i, m_indexBroken) = 0;
    }
  }

  data(i, m_indexBrokenNow) = 0;
}
//------------------------------------------------------------------------------
// void MohrCoulombMaxFractureWeighted::evaluateStepOne()
void MohrCoulombMaxFractureWeighted::evaluateStepTwo() {
  if (m_broken) {
    m_state = true;
  } else {
    m_state = false;
  }

  m_broken = false;
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFractureWeighted::exchangeBrokenParticlesMPI() {
#if USE_MPI
  vector<int> sendData;
  const vector<int> &neighbouringCores = m_grid->neighbouringCores();

  // Collecting all the send data
  for (auto broken_idList : m_brokenParticles) {
    const int id_i = broken_idList.first;
    const vector<int> &broken = broken_idList.second;

    sendData.push_back(id_i);
    sendData.push_back(broken.size());

    for (int b : broken)
      sendData.push_back(b);
  }
  const int myRank = MPI::COMM_WORLD.Get_rank();

  for (int toCore : neighbouringCores) {
    int nSendElements = sendData.size();
    int nRecieveElements;
    MPI_Status status;

    MPI_Sendrecv(&nSendElements, 1, MPI_INT, toCore, myRank, &nRecieveElements,
                 1, MPI_INT, toCore, toCore, MPI_COMM_WORLD, &status);

    int recieveData[nRecieveElements];

    MPI_Sendrecv(&sendData[0], nSendElements, MPI_INT, toCore, 1, &recieveData,
                 nRecieveElements, MPI_INT, toCore, 1, MPI_COMM_WORLD, &status);

    int i = 0;
    while (i < nRecieveElements) {
      const int id = recieveData[i++];
      int nBrokenBonds = recieveData[i++];
      for (int j = 0; j < nBrokenBonds; j++) {
        m_brokenParticles[id].push_back(recieveData[i++]);
      }
    }
  }
#endif
}
//------------------------------------------------------------------------------
}
