#include "adrmohrcoulombbondfracture.h"
#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
AdrMohrCoulombBondFracture::AdrMohrCoulombBondFracture(double mu, double C,
                                                       double T)
    : m_C(C), m_T(T), m_d(pow(sqrt(1 + mu * mu) + mu, 2)) {
  m_neededProperties = {pair<string, int>("stress", 1)};
  m_hasStepTwo = true;
}
//------------------------------------------------------------------------------
void AdrMohrCoulombBondFracture::initialize() {
  m_data = &m_particles->data();
  m_indexUnbreakable = m_particles->getParamId("unbreakable");
  m_indexConnected = m_particles->getPdParamId("connected");
  m_idToCol = &m_particles->idToCol();
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);

  m_state = false;
  m_maxPId = pair<int, int>(-1, -1);
  m_maxStress = std::numeric_limits<double>::min();

  switch (m_dim) {
  case 1:
    m_nStressElements = 1;
    m_ghostParameters = {"s_xx"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    break;
  case 2:
    m_nStressElements = 3;
    m_ghostParameters = {"s_xx", "s_yy", "s_xy"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    break;
  case 3:
    m_nStressElements = 6;
    m_ghostParameters = {"s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"};
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStress[3] = m_particles->registerParameter("s_zz");
    m_indexStress[4] = m_particles->registerParameter("s_xz");
    m_indexStress[5] = m_particles->registerParameter("s_yz");
    break;
  }
}
//------------------------------------------------------------------------------
void AdrMohrCoulombBondFracture::evaluateStepTwo(const int id_i, const int i) {
  // First calculating the total stress on an material point. The stress is
  // averaged
  // all its bonds, and the stress state on a bond is the mean of the
  // stress at each material point.
  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;
  mat &data = *m_data;
  const mat &R = m_particles->r();
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);
  double dr_ij[m_dim];

  if (m_dim == 2) {
    int counter = 0;
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol).at(id_j);

      if (data(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      double r_len = 0;
      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = R(j, d) - R(i, d);
        r_len += dr_ij[d] * dr_ij[d];
      }
      r_len = sqrt(r_len);

      double sx = 0.5 * (data(i, m_indexStress[0]) + data(j, m_indexStress[0]));
      double sy = 0.5 * (data(i, m_indexStress[1]) + data(j, m_indexStress[1]));
      double sxy =
          0.5 * (data(i, m_indexStress[2]) + data(j, m_indexStress[2]));
      //            double sx = max(fabs(data(i, m_indexStress[0])),
      //            fabs(data(j, m_indexStress[0])));
      //            double sy = max(fabs(data(i, m_indexStress[1])),
      //            fabs(data(j, m_indexStress[1])));
      //            double sxy = max(fabs(data(i, m_indexStress[2])),
      //            fabs(data(j, m_indexStress[2])));

      const double c = dr_ij[0] / r_len; // cos(theta)
      const double s = dr_ij[1] / r_len; // sin(theta)
      const double c2 = c * c;
      const double s2 = s * s;
      const double cs = c * s;
      const double s_n = sx * c2 + sy * s2 + 2 * sxy * cs;
      const double s_s = (sy - sx) * s * c + sxy * (c2 - s2);

      //            if(p_1 > m_T)
      //            {
      //                con.second[m_indexConnected] = 0;
      //                m_maxPId = pair<int, int>(id_i, counter);
      //            }

      //            if(m_d*p_1 - p_2 - m_C > 0)
      //            {
      //                con.second[m_indexConnected] = 0;
      //                m_maxPId = pair<int, int>(id_i, counter);
      ////                cout << "shearing/compression" << endl;
      //            }
      //            else
      //            if(p_1 > m_T)
      //            {
      //                con.second[m_indexConnected] = 0;
      //                m_maxPId = pair<int, int>(id_i, counter);
      //            }

      if (s_n >= m_T) {
        con.second[m_indexConnected] = 0;
        m_maxPId = pair<int, int>(id_i, counter);
        data(i, m_indexBrokenNow) = 1;
        /*
        const double theta_b = atan(dr_ij[1]/ dr_ij[0]);
        const double theta_p = 0.5*atan(2*sxy/(sx - sy));
        const double first = 0.5*(sx + sy);
        const double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);
        const double s1 = first + second;
        const double s2 = first - second;

        if(fabs(theta_b * theta_p)*180/M_PI > 30)
        {
        cout << id_i << " > " << id_j << " " << endl
             << s_n << "\t"<< s_s << "\t:" << theta_b*180/M_PI << endl
             << s1 << "\t" << 0.5*(s1 + s2) << "\t:" << theta_p*180/M_PI <<
        endl;
        }
        */
      }
      counter++;
    }
  }
}
//------------------------------------------------------------------------------
void AdrMohrCoulombBondFracture::evaluateStepTwo() {
  if (m_maxPId.first != -1) {
    m_state = true;
  } else {
    m_state = false;
  }

  m_maxPId = std::pair<int, int>(-1, -1);
  m_maxStress = std::numeric_limits<double>::min();
}
//------------------------------------------------------------------------------
}
