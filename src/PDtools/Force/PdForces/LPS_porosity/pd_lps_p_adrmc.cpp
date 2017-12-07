#include "pd_lps_p_adrmc.h"

#define USE_PRINCIPAL_STRESS 1
#define CALCULATE_NUMMERICAL_PRINCIPAL_STRESS 1

namespace PDtools {
//------------------------------------------------------------------------------
PD_LPS_porosity_adrmc::PD_LPS_porosity_adrmc(PD_Particles &particles, double m,
                                             double b, double phi, double C,
                                             double T, bool planeStress,
                                             bool analyticalM)
    : LPS_porosity_mc(particles, 0., m, b, phi, C, T, planeStress,
                      analyticalM) {
  particles.setNeedGhostVelocity(false);
  m_hasStaticModifier = true;
}
//------------------------------------------------------------------------------
void PD_LPS_porosity_adrmc::calculateForces(const int id, const int i) {
  const double theta_i = m_data(i, m_iTheta);
  const double m_i = m_data(i, m_iMass);
  const double a_i = m_data(i, m_iA);
  const double b_i = m_data(i, m_iA);

  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);

  const int nConnections = PDconnections.size();
  double dr0_ij[m_dim];
  double dr_ij[m_dim];
  _F.zeros();

  double thetaNew = 0;
  int nConnected = 0;

  //----------------------------------
  // TMP - standard stress calc from
  //    m_data(i, m_indexStress[0]) = 0;
  //    m_data(i, m_indexStress[1]) = 0;
  //    m_data(i, m_indexStress[2]) = 0;
  //----------------------------------

  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections[l_j];

    if (con.second[m_iConnected] <= 0.5)
      continue;

    const int id_j = con.first;
    const int j = m_idToCol_v[id_j];

    const double m_j = m_data(j, m_iMass);
    const double theta_j = m_data(j, m_iTheta);
    const double vol_j = m_data(j, m_iVolume);
    const double dr0 = con.second[m_iDr0];
    const double volumeScaling = con.second[m_iVolumeScaling];
    const double vol = vol_j * volumeScaling;
    const double w = weightFunction(dr0);
    const double a_j = m_data(j, m_iA);
    const double b_j = m_data(j, m_iB);

    double dr2 = 0;

    for (int d = 0; d < m_dim; d++) {
      dr0_ij[d] = m_r0(j, d) - m_r0(i, d);
      dr_ij[d] = m_r(j, d) - m_r(i, d);
      dr2 += dr_ij[d] * dr_ij[d];
    }

    const double dr = sqrt(dr2);
    const double ds = dr - dr0;
    double bond = (b_i * theta_i / m_i + b_j * theta_j / m_j) * dr0;
    bond += (a_i / m_i + a_j / m_j) * ds;
    bond *= w * vol / dr;
    thetaNew += w * dr0 * ds * vol;

    for (int d = 0; d < m_dim; d++) {
      m_F(i, d) += dr_ij[d] * bond;

      for (int d2 = 0; d2 < m_dim; d2++) {
        _F(d, d2) += w * dr_ij[d] * dr0_ij[d2] * vol;
      }
    }

    con.second[m_iStretch] = ds / dr0;
    //----------------------------------
    // TMP - standard stres calc from
    //        m_data(i, m_indexStress[0]) += 0.5*dr_ij[0]*dr_ij[0]*bond;
    //        m_data(i, m_indexStress[1]) += 0.5*dr_ij[1]*dr_ij[1]*bond;
    //        m_data(i, m_indexStress[2]) += 0.5*dr_ij[0]*dr_ij[1]*bond;
    //----------------------------------
    nConnected++;
  }

  if (nConnections <= 3) {
    m_data(i, m_iThetaNew) = 0;
  } else {
    m_data(i, m_iThetaNew) = m_dim / m_i * thetaNew;
  }

  //----------------------------------
  // TMP - standard stres calc from
  //----------------------------------
  //    if(nConnected > 5) {
  //        computeStress(id, i, nConnected);
  //    }
  computeStress(id, i, nConnected);
  //--------------------
  m_continueState = false;
}
//------------------------------------------------------------------------------
void PD_LPS_porosity_adrmc::evaluateStatic(int id, int i) {
  if (m_data(i, m_indexUnbreakable) >= 1)
    return;
#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
  arma::vec eigval(m_dim);
#endif
  vector<pair<int, vector<double>>> &PDconnections =
      m_particles.pdConnections(id);
  const double shearCrit = m_C0 - m_ks * m_T;
  //    const double s_crit = 3.*m_T/40.e9;

  bool broken = false;
  if (m_dim == 2) {
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      if (m_data(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      //            const double s = con.second[m_iStretch];

      const double sx = 0.5 * (m_data(i, m_indexStress[0]) + m_data(j, m_indexStress[0]));
      const double sy = 0.5 * (m_data(i, m_indexStress[1]) + m_data(j, m_indexStress[1]));
      const double sxy = 0.5 * (m_data(i, m_indexStress[2]) + m_data(j, m_indexStress[2]));

      const double first = 0.5 * (sx + sy);
      const double second = sqrt(0.25 * (sx - sy) * (sx - sy) + sxy * sxy);

      const double p_1 = first + second; // max
      const double p_2 = first - second; // min
#if USE_PRINCIPAL_STRESS
      const int criticalShear = p_2 <= m_C0 - m_ks * p_1;
#else
      const double shear_max = 0.5 * (p_1 - p_2);
      const double shear = shear_max * m_cos_theta;
      const double normal = 0.5 * (p_1 + p_2) + shear_max * m_sin_theta;
      const double criticalShear = shear > m_S0 - m_d * normal;
//            const double criticalShear = shear + m_d*normal - m_S0;
//            const double criticalTensile = p_1 - m_T;
#endif
      const int MC_valid = p_2 < shearCrit;
      const int criticalTensile = p_1 >= m_T;
      //            const int criticalTensile = s >= s_crit;

      if (MC_valid) {
        if (criticalShear) {
          m_data(i, m_indexBrokenNow) = 1;
          con.second[m_indexConnected] = 0;
          m_continueState = true;
          broken = true;
          //                    cout << "Shear\t " << id << " - " << id_j <<
          //                    "\tp1:" << p_1 << ", " << p_2 << endl;
        }
      }
      //            } else {
      //            }
      if (criticalTensile) {
        m_data(i, m_indexBrokenNow) = 1;
        con.second[m_indexConnected] = 0;
        m_continueState = true;
        broken = true;
        //                cout << "Tensile\t " << id << " - " << id_j <<
        //                "\tp1:" << p_1 << ", " << p_2 << endl;
      }
    }
  } else if (m_dim == 3) {
    arma::mat S(m_dim, m_dim);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int j = m_idToCol_v[id_j];

      if (m_data(j, m_indexUnbreakable) >= 1)
        continue;

      if (con.second[m_indexConnected] <= 0.5)
        continue;

      S(0, 0) = 0.5 * (m_data(i, m_indexStress[0]) + m_data(j, m_indexStress[0]));
      S(1, 1) = 0.5 * (m_data(i, m_indexStress[1]) + m_data(j, m_indexStress[1]));
      S(0, 1) = 0.5 * (m_data(i, m_indexStress[2]) + m_data(j, m_indexStress[2]));
      S(1, 0) = S(0, 1);
      S(2, 2) = 0.5 * (m_data(i, m_indexStress[3]) + m_data(j, m_indexStress[3]));
      S(0, 2) = 0.5 * (m_data(i, m_indexStress[4]) + m_data(j, m_indexStress[4]));
      S(2, 0) = S(0, 2);
      S(1, 2) = 0.5 * (m_data(i, m_indexStress[5]) + m_data(j, m_indexStress[5]));
      S(2, 1) = S(1, 2);

#if CALCULATE_NUMMERICAL_PRINCIPAL_STRESS
      arma::eig_sym(eigval, S);
      const double p_1 = eigval(2);
      const double p_2 = eigval(0);
#else
      const double I1 = S(0, 0) + S(1, 1) + S(2, 2);
      const double I2 = S(0, 0) * S(1, 1) + S(1, 1) * S(2, 2) +
                        S(3, 3) * S(0, 0) - pow(S(0, 1), 2) - pow(S(1, 2), 2) -
                        pow(S(0, 2), 2);
      const double I3 = S(0, 0) * S(1, 1) * S(2, 2) -
                        S(0, 0) * pow(S(1, 2), 2) - S(1, 1) * pow(S(0, 2), 2) -
                        S(2, 2) * pow(S(0, 1), 2) +
                        2 * S(0, 1) * S(1, 2) * S(0, 2);
      const double phi =
          1. / 3. * acos(0.5 * (2 * pow(I1, 3) - 9 * I1 * I2 + 27 * I3) /
                         pow(pow(I1, 2) - 3 * I2, 1.5));

      const double core = 2. / 3. * (sqrt(I1 * I1 - 3 * I2));
      const double s1 = I1 / 3. + core * cos(phi);
      const double s2 = I1 / 3. + core * cos(phi - 2. * M_PI / 3.);
      const double s3 = I1 / 3. + core * cos(phi - 4. * M_PI / 3.);

      double p_1 = s1;
      double p_2 = s2;

      if (s2 > p_1) {
        p_1 = s2;
        p_2 = s1;
      }
      if (s3 > p_1) {
        p_1 = s3;
      }

      if (p_2 < s3) {
        p_2 = s3;
      }
#endif
      const int criticalShear = p_2 <= m_C0 - m_ks * p_1;
      const int MC_valid = p_2 < shearCrit;
      const int criticalTensile = p_1 >= m_T;

      if (MC_valid) {
        if (criticalShear) {
          m_data(i, m_indexBrokenNow) = 1;
          con.second[m_indexConnected] = 0;
          broken = true;
        }
      } else {
        if (criticalTensile) {
          m_data(i, m_indexBrokenNow) = 1;
          con.second[m_indexConnected] = 0;
          broken = true;
        }
      }
    }
  }

  if (broken) {
    updateWeightedVolume(id, i);
    computeK(id, i);
    m_data(i, m_indexBrokenNow) = 0;
  }
}

//------------------------------------------------------------------------------
void PD_LPS_porosity_adrmc::evaluateStepTwo(int id, int i) {
  (void)id;
  (void)i;
}
//------------------------------------------------------------------------------
}
