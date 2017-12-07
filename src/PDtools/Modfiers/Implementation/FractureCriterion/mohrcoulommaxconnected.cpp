#include "mohrcoulommaxconnected.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
MohrCoulomMaxConnected::MohrCoulomMaxConnected(double mu, double C, double T)
    : m_C(C), m_T(T) {
  m_phi = mu * M_PI / 180.;

  m_d = tan(m_phi);
  m_neededProperties = {pair<string, int>("stress", 1)};

  m_hasStepOne = true;
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void MohrCoulomMaxConnected::registerParticleParameters() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexCompute = m_particles->registerPdParameter("compute");
  m_indexBroken = m_particles->registerParameter("broken", 0);
  m_indexBrokenId = m_particles->registerParameter("brokendId", 0);
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_idToCol = &m_particles->getIdToCol_v();

  switch (m_dim) {
  case 1:
    m_ghostParameters = {"nx"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexNormal[0] = m_particles->registerParameter("nx");
    break;
  case 2:
    m_ghostParameters = {"nx", "ny"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexNormal[0] = m_particles->registerParameter("nx");
    m_indexNormal[1] = m_particles->registerParameter("ny");
    break;
  case 3:
    m_ghostParameters = {"nx", "ny", "nz"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStress[3] = m_particles->registerParameter("s_zz");
    m_indexStress[4] = m_particles->registerParameter("s_xz");
    m_indexStress[5] = m_particles->registerParameter("s_yz");
    m_indexNormal[0] = m_particles->registerParameter("nx");
    m_indexNormal[1] = m_particles->registerParameter("ny");
    m_indexNormal[2] = m_particles->registerParameter("nz");
    break;
  }
  m_ghostParameters.push_back("brokendId");
}
//------------------------------------------------------------------------------
void MohrCoulomMaxConnected::initialize() {
  m_broken = false;
  m_cosTheta = cos(M_PI / 2. + m_phi);
  m_sinTheta = sin(M_PI / 2. + m_phi);
}
//------------------------------------------------------------------------------
void MohrCoulomMaxConnected::evaluateStepOne(const int id_i, const int i) {
  mat &data = *m_data;
  const mat &r = m_particles->r();

  if (data(i, m_indexUnbreakable) >= 1)
    return;

  const int broken_i = data(i, m_indexBroken);
  double n_i[m_dim];
  double n_j[m_dim];
  double dr_ij[m_dim];

  if (broken_i) {
    for (int d = 0; d < m_dim; d++) {
      n_i[d] = data(i, m_indexNormal[d]);
    }
  }

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);

  for (auto &con : PDconnections) {
    const int id_j = con.first;
    const int j = (*m_idToCol)[id_j];

    if (data(j, m_indexUnbreakable) >= 1)
      continue;

    if (con.second[m_indexConnected] <= 0.5)
      continue;

    if (broken_i) {

      double dotproduct = 0;
      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = r(j, d) - r(i, d);
        dotproduct += n_i[d] * dr_ij[d];
      }

      if (dotproduct > 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }
      const int bId = data(i, m_indexBrokenId);
      if (bId == id_j) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }

      continue;
    }

    const int broken_j = data(j, m_indexBroken);
    if (broken_j > 0) {
      const int bId = data(i, m_indexBrokenId);

      if (bId == id_i) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }

      double dotproduct = 0;
      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = r(i, d) - r(j, d);
        n_j[d] = data(j, m_indexNormal[d]);
        dotproduct += n_j[d] * dr_ij[d];
      }

      if (dotproduct > 0) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;
      }

      //--------------------------------------------------------------
      // Cheking all the othder bonds
      //--------------------------------------------------------------
      n_j[0] = -data(j, m_indexNormal[1]);
      n_j[1] = data(j, m_indexNormal[0]);
      double x1, x2, x3, x4;
      double y1, y2, y3, y4;
      x1 = r(j, 0);
      y1 = r(j, 1);
      x2 = r(j, 0) + n_j[0];
      y2 = r(j, 1) + n_j[1];
      x3 = r(i, 0);
      y3 = r(i, 1);

      const double center_x = r(j, 0) + 0.5 * n_j[0];
      const double center_y = r(j, 1) + 0.5 * n_j[1];
      const double r2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

      vector<pair<int, vector<double>>> &PDconnections2 =
          m_particles->pdConnections(id_i);
      for (auto &con2 : PDconnections2) {
        if (con2.first == id_j) {
          continue;
        }
        const int id_k = con2.first;
        const int k = (*m_idToCol).at(id_k);
        x4 = r(k, 0);
        y4 = r(k, 1);

        // If the denominator is close to zero the two lines are parallell of
        // coincident
        const double threshold = 1e-7;
        const double det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (fabs(det) > threshold)
          continue;

        // Finding the intersection between the two lines
        const double px = ((x1 * y2 - y1 * x2) * (x3 - x4) -
                           (x1 - x2) * (x3 * y4 - y3 * x4)) /
                          det;
        const double py = ((x1 * y2 - y1 * x2) * (y3 - y4) -
                           (y1 - y2) * (x3 * y4 - y3 * x4)) /
                          det;
        const double dist = (px - center_x) * (px - center_x) +
                            (py - center_y) * (py - center_y);

        if (dist < r2) {
          double dotproduct = (px - r(i, 0)) * (px - r(k, 0));
          dotproduct += (py - r(i, 1)) * (py - r(k, 1));
          if (dotproduct < 0) {
            con2.second[m_indexConnected] = 0;
            vector<pair<int, vector<double>>> &PDconnections_k =
                m_particles->pdConnections(id_k);
            for (auto &con_k : PDconnections_k) {
              if (con_k.first == id_i) {
                con_k.second[m_indexConnected] = 0;
                data(i, m_indexBrokenNow) = 1;
              }
            }
          }
        }
      }
    }
  }
}
//------------------------------------------------------------------------------
void MohrCoulomMaxConnected::evaluateStepTwo(const int id_i, const int i) {
  (void)id_i;
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  mat &data = *m_data;
  const mat &r = m_particles->r();

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

    int broken = 0;

    if (fabs(shear) >= fabs(m_C - m_d * normal) && normal <= 0) {
      data(i, m_indexBroken) = 2;
      m_broken = true;
      broken = 2;
    } else if (p_2 >= m_T && normal > 0) {
      data(i, m_indexBroken) = 1;
      m_broken = true;
      broken = 1;
    } else {
      data(i, m_indexBroken) = 0;
    }

    double closest = std::numeric_limits<double>::max();
    if (broken > 0) {
      vector<pair<int, vector<double>>> &PDconnections =
          m_particles->pdConnections(id_i);
      int indexMax = -1;
      int colMax = -1;

      for (auto &con : PDconnections) {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if (data(j, m_indexUnbreakable) >= 1)
          continue;

        if (con.second[m_indexConnected] <= 0.5)
          continue;

        sx = data(j, m_indexStress[0]);
        sy = data(j, m_indexStress[1]);
        sxy = data(j, m_indexStress[2]);

        const double first_j = 0.5 * (sx + sy);
        const double second_j = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

        const double pj_2 = first_j + second_j; // max
        const double pj_1 = first_j - second_j; // min

        //                const double shear_j = 0.5*(pj_1 - pj_2)*m_sinTheta;
        //                const double normal_j = 0.5*(pj_1 + pj_2) + 0.5*(pj_1
        //                - pj_2)*m_cosTheta;

        if (broken == 1) {
          const double dp = p_2 - pj_2;
          if (dp < closest) {
            closest = dp;
            indexMax = id_j;
            colMax = j;
          }
        }
      }

      //            if(indexMax == -1)
      //                continue;
      data(i, m_indexBrokenId) = indexMax;
      //            data(i, m_indexNormal[0]) = -(r(i, 1) - r(colMax, 1));
      //            data(i, m_indexNormal[1]) = r(i, 0) - r(colMax, 0);
      data(i, m_indexNormal[0]) = r(i, 0) - r(colMax, 0);
      data(i, m_indexNormal[1]) = -(r(i, 1) - r(colMax, 1));
    }
  }
}
//------------------------------------------------------------------------------
}
