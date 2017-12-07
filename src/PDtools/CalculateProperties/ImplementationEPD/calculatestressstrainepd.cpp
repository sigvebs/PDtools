#include "calculatestressstrainepd.h"

#include "PDtools/Elements/pd_element.h"
#include "PDtools/Force/force.h"
#include "Particles/pd_particles.h"

namespace PDtools {
//------------------------------------------------------------------------------
CalculateStressStrainEPD::CalculateStressStrainEPD(vector<Force *> &forces,
                                                   double E, double nu,
                                                   double delta,
                                                   bool planeStress,
                                                   bool greenLagrange)
    : CalculateProperty("stress"), m_E(E), m_nu(nu), m_delta(delta),
      m_planeStress(planeStress) {
  for (Force *force : forces) {
    m_forces.push_back(force);
  }

  m_greenStrain = greenLagrange;
  if (m_greenStrain) {
    m_smallStrain = false;
  } else {
    m_smallStrain = true;
  }
}
//------------------------------------------------------------------------------
void CalculateStressStrainEPD::initialize() {
  m_particles->setNeedGhostR0(true);
  m_indexConnected = m_particles->getPdParamId("connected");
  //    m_indexVolume = m_particles->getParamId("volume");
  //    m_indexVolumeScaling = m_particles->getPdParamId("volumeScaling");
  //    m_indexDr0 = m_particles->getPdParamId("dr0");
  m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);
  m_iOverlap = m_particles->getPdParamId("overlap");

  switch (m_dim) {
  case 1:
    m_nStressStrainElements = 1;
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    m_indexK[0] = m_particles->registerParameter("K_x");
    break;
  case 2:
    m_nStressStrainElements = 3;
    m_indexK[0] = m_particles->registerParameter("K_x");
    m_indexK[1] = m_particles->registerParameter("K_y");
    m_indexK[2] = m_particles->registerParameter("K_xy");
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    m_indexStrain[1] = m_particles->registerParameter("e_yy");
    m_indexStrain[2] = m_particles->registerParameter("e_xy");
    break;
  case 3:
    m_nStressStrainElements = 6;
    m_indexK[0] = m_particles->registerParameter("K_x");
    m_indexK[1] = m_particles->registerParameter("K_y");
    m_indexK[2] = m_particles->registerParameter("K_xy");
    m_indexK[3] = m_particles->registerParameter("K_zz");
    m_indexK[4] = m_particles->registerParameter("K_xz");
    m_indexK[5] = m_particles->registerParameter("K_yz");
    m_indexStress[0] = m_particles->registerParameter("s_xx");
    m_indexStress[1] = m_particles->registerParameter("s_yy");
    m_indexStress[2] = m_particles->registerParameter("s_xy");
    m_indexStress[3] = m_particles->registerParameter("s_zz");
    m_indexStress[4] = m_particles->registerParameter("s_xz");
    m_indexStress[5] = m_particles->registerParameter("s_yz");
    m_indexStrain[0] = m_particles->registerParameter("e_xx");
    m_indexStrain[1] = m_particles->registerParameter("e_yy");
    m_indexStrain[2] = m_particles->registerParameter("e_xy");
    m_indexStrain[3] = m_particles->registerParameter("e_zz");
    m_indexStrain[4] = m_particles->registerParameter("e_xz");
    m_indexStrain[5] = m_particles->registerParameter("e_yz");
    break;
  }

  // Computing the shape tensor
  const ivec &colToId = m_particles->colToId();
  const int nParticles = m_particles->nParticles();

  for (int i = 0; i < nParticles; i++) {
    const int id_i = colToId(i);
    computeK(id_i, i);
  }
}
//------------------------------------------------------------------------------
void CalculateStressStrainEPD::clean() {}
//------------------------------------------------------------------------------
void CalculateStressStrainEPD::update() {
  const ivec &colToId = m_particles->colToId();
//  const auto &idToCol = m_particles->idToCol();
  const int nParticles = m_particles->nParticles();
  const mat &r = m_particles->r();
  const mat &r0 = m_particles->r0();
  mat &data = m_particles->data();
  vector<PD_quadElement> &quadElements = m_particles->getQuadElements();
  unordered_map<int, int> &idToElement = m_particles->getIdToElement();

  mat F = zeros(m_dim, m_dim);
  //    mat E = zeros(m_dim, m_dim);
  mat K = zeros(m_dim, m_dim);
  mat P = zeros(m_dim, m_dim);
  mat strain = zeros(m_dim, m_dim);

  double dr_ij[m_dim];
  double dr0_ij[m_dim];

  for (int i = 0; i < nParticles; i++) {
    const int id_i = colToId(i);
    if (data(i, m_indexBrokenNow)) {
      computeK(id_i, i);
      data(i, m_indexBrokenNow) = 0;
    }

    if (m_dim >= 2) {
      K(0, 0) = data(i, m_indexK[0]);
      K(1, 1) = data(i, m_indexK[1]);
      K(0, 1) = data(i, m_indexK[2]);
      K(1, 0) = K(0, 1);
    }

    F.zeros();
    vector<pair<int, vector<double>>> &PDconnections_i =
        m_particles->pdConnections(id_i);
    const int nConnections = PDconnections_i.size();
    int nConnected = 0;

    for (int l_j = 0; l_j < nConnections; l_j++) {
      const auto &con_i = PDconnections_i[l_j];
      if (con_i.second[m_indexConnected] <= 0.5)
        continue;

      const int polygon_id = con_i.first;
      const int polygon_i = idToElement.at(polygon_id);

      PD_quadElement &element = quadElements[polygon_i];
      const mat &gaussPoints_initial =
          element.guassianQuadraturePoints_initial();
      const mat &gaussPoints = element.guassianQuadraturePoints();
      const vec &gaussWeights = element.guassianQuadratureWeights();

      const int nIntegrationPoints = gaussPoints.n_rows;

      for (int j = 0; j < nIntegrationPoints; j++) {
        const double w_j = gaussWeights[j];
        double dr0 = 0;
        for (int d = 0; d < m_dim; d++) {
          dr_ij[d] = gaussPoints(j, d) - r(i, d);
          dr0_ij[d] = gaussPoints_initial(j, d) - r0(i, d);
          dr0 += dr0_ij[d] * dr0_ij[d];
        }

#if OVERLAP
        const double overlap = con_i.second[m_iOverlap];
        dr0 = sqrt(dr0);
        const double w = overlap * weightFunction(dr0);
#else
        if (dr0 > m_delta * m_delta)
          continue;
        dr0 = sqrt(dr0);
        const double w = weightFunction(dr0);
#endif

        for (int d = 0; d < m_dim; d++) {
          for (int d2 = 0; d2 < m_dim; d2++) {
            F(d, d2) += w * dr_ij[d] * dr0_ij[d2] * w_j;
          }
        }
        nConnected++;
      }
    }

    if (nConnected <= 3) {
      data(i, m_indexStrain[0]) = 0;
      data(i, m_indexStrain[1]) = 0;
      data(i, m_indexStrain[2]) = 0;
      data(i, m_indexStress[0]) = 0;
      data(i, m_indexStress[1]) = 0;
      data(i, m_indexStress[2]) = 0;
      continue;
    } else {
      F = F * K; // K = inv(K);
      if (m_smallStrain) {
        strain = 0.5 * (F.t() + F);
        strain(0, 0) -= 1.;
        strain(1, 1) -= 1.;
      } else {
        strain = 0.5 * F.t() * F;
        strain(0, 0) -= 0.5;
        strain(1, 1) -= 0.5;
      }
    }

    data(i, m_indexStrain[0]) = strain(0, 0);
    data(i, m_indexStrain[1]) = strain(1, 1);
    data(i, m_indexStrain[2]) = strain(0, 1);

    // Assuming linear elasticity
    if (m_dim == 2) {
      // Constituent model, linear elastic
      // Computing the second PK stress
      if (m_planeStress) {
        const double a = m_E / (1. - m_nu * m_nu);
        P(0, 0) = a * (strain(0, 0) + m_nu * strain(1, 1));
        P(1, 1) = a * (strain(1, 1) + m_nu * strain(0, 0));
        P(0, 1) = a * (1 - m_nu) * strain(0, 1);
        P(1, 0) = P(0, 1);
      } else // Plane strain
      {
        const double a = m_E / ((1. + m_nu) * (1. - 2 * m_nu));
        P(0, 0) = a * ((1. - m_nu) * strain(0, 0) + m_nu * strain(1, 1));
        P(1, 1) = a * ((1. - m_nu) * strain(1, 1) + m_nu * strain(0, 0));
        P(0, 1) = a * 0.5 * (1. - 2. * m_nu) * strain(0, 1);
        P(1, 0) = P(0, 1);
      }

      // Converting PK stress to Cauchy stress
      if (m_greenStrain) {
        const double detF = 1. / det(F);
        P = detF * F * P * F.t();
      }
      data(i, m_indexStress[0]) = P(0, 0);
      data(i, m_indexStress[1]) = P(1, 1);
      data(i, m_indexStress[2]) = P(0, 1);
    }

    if (m_dim == 3) {
      //            F(2,2) -= 0.5;
      //            data(i, m_indexStrain[3]) = F(2,2);
      //            data(i, m_indexStrain[4]) = F(0,2);
      //            data(i, m_indexStrain[5]) = F(1,2);
    }
  }
}
//------------------------------------------------------------------------------
void CalculateStressStrainEPD::computeK(int id, int i) {
  const mat &r0 = m_particles->r0();
  mat &data = m_particles->data();
  mat K = zeros(m_dim, m_dim);
  vector<pair<int, vector<double>>> &PDconnections_i =
      m_particles->pdConnections(id);
  const int nConnections = PDconnections_i.size();
  double dr0_ij[m_dim];

  vector<PD_quadElement> &quadElements = m_particles->getQuadElements();
  unordered_map<int, int> &idToElement = m_particles->getIdToElement();

  int nConnected = 0;
  for (int l_j = 0; l_j < nConnections; l_j++) {
    auto &con = PDconnections_i[l_j];
    if (con.second[m_indexConnected] <= 0.5)
      continue;

    const int polygon_id = con.first;
//    const double overlap = con.second[m_iOverlap];
    const int polygon_i = idToElement.at(polygon_id);

    PD_quadElement &element = quadElements[polygon_i];
    const mat &gaussPoints_initial = element.guassianQuadraturePoints_initial();
    const mat &gaussPoints = element.guassianQuadraturePoints();
    const vec &gaussWeights = element.guassianQuadratureWeights();

    const int nIntegrationPoints = gaussPoints.n_rows;

    for (int j = 0; j < nIntegrationPoints; j++) {
      const double w_j = gaussWeights[j];

      double dr2_0 = 0;
      for (int d = 0; d < m_dim; d++) {
        dr0_ij[d] = gaussPoints_initial(j, d) - r0(i, d);
        dr2_0 += dr0_ij[d] * dr0_ij[d];
      }
#if OVERLAP
      dr0 = sqrt(dr0);
      const double w = overlap * weightFunction(dr0);
#else
      if (dr2_0 > m_delta * m_delta)
        continue;
      dr2_0 = sqrt(dr2_0);
      const double w = weightFunction(dr2_0);
#endif

      for (int d = 0; d < m_dim; d++) {
        for (int d2 = 0; d2 < m_dim; d2++) {
          K(d, d2) += w * dr0_ij[d] * dr0_ij[d2] * w_j;
        }
      }
      nConnected++;
    }
  }

  if (nConnected <= 3) {
    if (m_dim >= 2) {
      data(i, m_indexK[0]) = 0;
      data(i, m_indexK[1]) = 0;
      data(i, m_indexK[2]) = 0;
    }
    if (m_dim == 3) {
      data(i, m_indexK[3]) = 0;
      data(i, m_indexK[4]) = 0;
      data(i, m_indexK[5]) = 0;
    }
  } else {
    //        K = inv_sympd(K);
    K = inv(K);
  }

  if (m_dim >= 2) {
    data(i, m_indexK[0]) = K(0, 0);
    data(i, m_indexK[1]) = K(1, 1);
    data(i, m_indexK[2]) = K(0, 1);
  }
  if (m_dim == 3) {
    data(i, m_indexK[3]) = K(2, 2);
    data(i, m_indexK[4]) = K(0, 2);
    data(i, m_indexK[5]) = K(1, 2);
  }
}
//------------------------------------------------------------------------------
}
