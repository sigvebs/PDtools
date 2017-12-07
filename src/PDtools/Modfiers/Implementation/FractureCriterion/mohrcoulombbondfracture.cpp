#include "mohrcoulombbondfracture.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
MohrCoulombBondFracture::MohrCoulombBondFracture(double mu, double C, double T)
    : m_C(C), m_T(T) {
  m_d = tan(mu * M_PI / 180.);
  m_neededProperties = {pair<string, int>("stress", 1)};
  m_hasStepOne = true;
}
//------------------------------------------------------------------------------
MohrCoulombBondFracture::~MohrCoulombBondFracture() {}
//------------------------------------------------------------------------------
void MohrCoulombBondFracture::registerParticleParameters() {
  m_data = &m_particles->data();
  m_r = &m_particles->r();
  m_indexUnbreakable = m_particles->registerParameter("unbreakable");
  m_indexConnected = m_particles->registerPdParameter("connected");
  m_indexCompute = m_particles->registerPdParameter("compute");
  m_idToCol = &m_particles->getIdToCol_v();
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);

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
    cerr << "MohrCoulombBondFracture: not implemented in 3D." << endl;
    exit(1);
    break;
  }
}
//------------------------------------------------------------------------------
void MohrCoulombBondFracture::initialize() {}
//------------------------------------------------------------------------------
void MohrCoulombBondFracture::evaluateStepOne(const int id_i, const int i) {
  // First calculating the total stress on an material point. The stress is
  // averaged
  // all its bonds, and the stress state on a bond is the mean of the
  // stress at each material point.

  if ((*m_data)(i, m_indexUnbreakable) >= 1)
    return;

  const mat &R = m_particles->r();
  mat &data = *m_data;
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles->pdConnections(id_i);
  //    arma::mat S_i(m_dim, m_dim);
  //    arma::mat S(m_dim, m_dim);
  double dr_ij[m_dim];

  if (m_dim == 2) {
    //        S_i(0, 0) = data(i, m_indexStress[0]);
    //        S_i(1, 1) = data(i, m_indexStress[1]);
    //        S_i(0, 1) = data(i, m_indexStress[2]);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = (*m_idToCol)[id_j];

      if ((*m_data)(j, m_indexUnbreakable) >= 1)
        continue;
      if (con.second[m_indexConnected] <= 0.5)
        continue;

      // Finding the bond angle
      double r_len = 0;
      for (int d = 0; d < m_dim; d++) {
        dr_ij[d] = R(j, d) - R(i, d);
        r_len += dr_ij[d] * dr_ij[d];
      }
      r_len = sqrt(r_len);
      //            const double theta = atan2(dr_ij[1], dr_ij[0]);

      double sx = 0.5 * (data(i, m_indexStress[0]) + data(j, m_indexStress[0]));
      double sy = 0.5 * (data(i, m_indexStress[1]) + data(j, m_indexStress[1]));
      double sxy = 0.5 * (data(i, m_indexStress[2]) + data(j, m_indexStress[2]));
      //            double sx = max(fabs(data(i, m_indexStress[0])),
      //            fabs(data(j, m_indexStress[0])));
      //            double sy = max(fabs(data(i, m_indexStress[1])),
      //            fabs(data(j, m_indexStress[1])));
      //            double sxy = max(fabs(data(i, m_indexStress[2])),
      //            fabs(data(j, m_indexStress[2])));
      //            const double c = cos(theta);
      //            const double s = sin(theta);
      const double c = dr_ij[0] / r_len; // cos(theta)
      const double s = dr_ij[1] / r_len; // sin(theta)
      const double c2 = c * c;
      const double s2 = s * s;
      const double cs = c * s;
      const double s_n = sx * c2 + sy * s2 + 2 * sxy * cs;
      const double s_s = (sy - sx) * s * c + sxy * (c2 - s2);

      if (s_n >= m_T) {
        con.second[m_indexConnected] = 0;
        data(i, m_indexBrokenNow) = 1;

        const double theta_b = atan(dr_ij[1] / dr_ij[0]);
        const double theta_p = 0.5 * atan(2 * sxy / (sx - sy));
        const double first = 0.5 * (sx + sy);
        const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);
        const double s1 = first + second;
        const double s2 = first - second;

        if (fabs(theta_b * theta_p) * 180 / M_PI > 30) {
          cout << id_i << " > " << id_j << " " << endl
               << s_n << "\t" << s_s << "\t:" << theta_b * 180 / M_PI << endl
               << s1 << "\t" << 0.5 * (s1 + s2) << "\t:" << theta_p * 180 / M_PI
               << endl;
        }
      } else if (s_n < 0) {
        if (fabs(s_s) >= m_C - m_d * s_n) {
          con.second[m_indexConnected] = 0;
          data(i, m_indexBrokenNow) = 1;
        }
      }

      //            sy =  sx*s2 + sy*c2 - 2*sxy*cs;

      //            if(sxy > fabs(-m_C + m_d*s_compressive))
      //            {
      //                con.second[m_indexConnected] = 0;
      //                cout << id_i << " -> (" << id_j << ") " << sx << " "<<
      //                sy << " t:" << theta_i*180/M_PI << endl;
      //                cout << "sx_i:" << data(i, m_indexStress[0]) << " sx_j:"
      //                << data(j, m_indexStress[0]) << endl;
      //                cout << "sy_i:" << data(i, m_indexStress[1]) << " sy_j:"
      //                << data(j, m_indexStress[1]) << endl;
      //            }
    }

  } else if (m_dim == 3) {
  }
}
//------------------------------------------------------------------------------
}
