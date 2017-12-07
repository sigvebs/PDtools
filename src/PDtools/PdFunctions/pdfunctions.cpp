#include "pdfunctions.h"

#include "Grid/grid.h"
#include "PDtools/Force/force.h"
#include "Particles/pd_particles.h"
#include <boost/algorithm/string.hpp>
#include <limits>
#include <math.h>
#include <unordered_map>

#ifdef USE_MPI
#include <mpi.h>
#endif

#define USE_EXTENDED_RANGE_RADIUS 0
#define USE_EXTENDED_RANGE_LC 1
//#define USE_EXTENDED_RANGE_LC 0

namespace PDtools {
//------------------------------------------------------------------------------
void setPdConnections(PD_Particles &particles, Grid &grid, double delta,
                      double lc) {
#if USE_EXTENDED_RANGE_LC == 0
  (void)lc;
#endif
  const unordered_map<int, GridPoint> &gridpoints = grid.gridpoints();
  const vector<int> &mygridPoints = grid.myGridPoints();
  const mat &R = particles.r();
  const mat &data = particles.data();
#if USE_EXTENDED_RANGE_RADIUS
  const int indexRadius = particles.getParamId("radius");
#endif
  // The order is important!
  particles.registerPdParameter("dr0");
  particles.registerPdParameter("connected");
  const int iGroupId = particles.getParamId("groupId");

//    int nBonds = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < mygridPoints.size(); i++) {
    double dx, dy, dz;
    int gridId = mygridPoints.at(i);
    const GridPoint &gridPoint = gridpoints.at(gridId);

    for (const pair<int, int> &idCol_i : gridPoint.particles()) {
      int id_i = idCol_i.first;
      int col_i = idCol_i.second;
      const vec &r_i = R.row(col_i).t();
      const int group_i = data(col_i, iGroupId);
      unordered_map<int, vector<double>> connections;
      vector<pair<int, vector<double>>> connectionsVector;

#if USE_EXTENDED_RANGE_RADIUS
      const double radius_i = data(col_i, indexRadius);
#elif USE_EXTENDED_RANGE_LC
      const double radius_i = 0.5 * lc;
#else
      const double radius_i = 0.;
#endif
      for (const pair<int, int> &idCol_j : gridPoint.particles()) {
        const int id_j = idCol_j.first;
        const int col_j = idCol_j.second;
        if (id_i == id_j)
          continue;

        const int group_j = data(col_j, iGroupId);
        if (group_i != group_j)
          continue;

#if USE_EXTENDED_RANGE_RADIUS
        const double radius_j = data(col_j, indexRadius);
#elif USE_EXTENDED_RANGE_LC
        const double radius_j = 0.5 * lc;
#else
        const double radius_j = 0.;
#endif

        const vec &r_j = R.row(col_j).t();
        dx = r_i(0) - r_j(0);
        dy = r_i(1) - r_j(1);
        dz = r_i(2) - r_j(2);

        const double drSquared = dx * dx + dy * dy + dz * dz;
        const double dr = sqrt(drSquared);

        if (delta >= dr - radius_i && delta >= dr - radius_j) {
          vector<double> connectionData;

          connectionData.push_back(dr);
          connectionData.push_back(1.0); // Connected
          connections[id_j] = connectionData;
          connectionsVector.push_back(
              pair<int, vector<double>>(id_j, connectionData));
        }
      }

      // Neighbouring cells
      const vector<GridPoint *> &neighbours = gridPoint.neighbours();

      for (const GridPoint *neighbour : neighbours) {
        for (const pair<int, int> &idCol_j : neighbour->particles()) {
          const int col_j = idCol_j.second;
          const int group_j = data(col_j, iGroupId);
          if (group_i != group_j)
            continue;

#if USE_EXTENDED_RANGE_RADIUS
          const double radius_j = data(col_j, indexRadius);
          const double l_delta = delta + 0.5 * (radius_i + radius_j);
#elif USE_EXTENDED_RANGE_LC
          const double radius_j = 0.5 * lc;
#else
          const double radius_j = 0.;
#endif

          const vec &r_j = R.row(col_j).t();
          dx = r_i(0) - r_j(0);
          dy = r_i(1) - r_j(1);
          dz = r_i(2) - r_j(2);

          const double drSquared = dx * dx + dy * dy + dz * dz;
          const double dr = sqrt(drSquared);

          if (delta >= dr - radius_i && delta >= dr - radius_j) {
            int id_j = idCol_j.first;
            vector<double> connectionData;
            connectionData.push_back(dr);
            connectionData.push_back(1.0); // Connected
            connections[id_j] = connectionData;
            connectionsVector.push_back(
                pair<int, vector<double>>(id_j, connectionData));
          }
        }
      }

#ifdef USE_OPENMP
#pragma omp critical
      { particles.setPdConnections(id_i, connectionsVector); }
#else
      particles.setPdConnections(id_i, connectionsVector);
#endif
    }
  }
}
//------------------------------------------------------------------------------
void applyVolumeCorrection(PD_Particles &particles, double delta, double lc,
                           int dim) {
#if USE_EXTENDED_RANGE_LC == 0
  (void)lc;
#endif
  const mat &data = particles.data();
  const ivec &idToCol = particles.getIdToCol_v();
  const ivec &colToId = particles.colToId();
  const int indexDr0 = particles.getPdParamId("dr0");
  //    const int indexRadius = particles.getParamId("radius");
  const int indexVolumeScaling = particles.getPdParamId("volumeScaling");
  const int indexVolume = particles.getParamId("volume");

  // http://mathworld.wolfram.com/Circle-CircleIntersection.html
  if (dim == 2) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const int id_i = colToId(i);
#if USE_EXTENDED_RANGE_RADIUS
      const double r_i = data(i, indexRadius);
#elif USE_EXTENDED_RANGE_LC
      const double r_i = 0.5 * lc;
#else
      const double r_i = 0.;
#endif

      vector<pair<int, vector<double>>> &PDconnections =
          particles.pdConnections(id_i);
      double vol_delta = 0;

      for (auto &con : PDconnections) {
        const double dr = con.second[indexDr0];
#if USE_EXTENDED_RANGE_RADIUS
        const int id_j = con.first;
        const int col_j = idToCol[id_j];
        const double r_j = data(col_j, indexRadius);
#elif USE_EXTENDED_RANGE_LC
        const double r_j = 0.5 * lc;
#else
        const double r_j = 0.;
#endif
        const double rc_j = delta - r_j;
        double v1 = 1.;

        if (dr > rc_j) {
          const double d = dr;
          const double R = delta;
          const double r = r_j;

          const double d1 = 0.5 * (d * d - r * r + R * R) / d;
          const double d2 = 0.5 * (d * d + r * r - R * R) / d;
          const double A1 = R * R * acos(d1 / R) - d1 * sqrt(R * R - d1 * d1);
          const double A2 = r * r * acos(d2 / r) - d2 * sqrt(r * r - d2 * d2);

          const double A = A1 + A2;
          const double A_a = M_PI * r * r;

          v1 = A / A_a;
        }

        double volumeCorrection = 1.0;
        volumeCorrection = v1;

        con.second[indexVolumeScaling] = volumeCorrection;
        const double vol_i = data(i, indexVolume);
        vol_delta += vol_i * volumeCorrection;
      }
    }
  } else if (dim == 3) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const int pId = colToId(i);
      vector<pair<int, vector<double>>> &PDconnections =
          particles.pdConnections(pId);
      double vol_delta = 0;

      for (auto &con : PDconnections) {
        const double dr = con.second[indexDr0];
#if USE_EXTENDED_RANGE_RADIUS
        const int id_j = con.first;
        const int col_j = idToCol[id_j];
        const double radius_j = data(col_j, indexRadius);
#elif USE_EXTENDED_RANGE_LC
        const double radius_j = 0.5 * lc;
#else
        const double radius_j = data(col_j, indexRadius);
#endif
        const double rc = delta - radius_j;

        double volumeCorrection = 1.0;
        if (dr > rc) {
          volumeCorrection = 0.5 * (delta + radius_j - dr) / radius_j;
        }
        con.second[indexVolumeScaling] = volumeCorrection;
        const double vol_i = data(i, indexVolume);
        vol_delta += vol_i * volumeCorrection;
      }
    }
  }
}
//------------------------------------------------------------------------------
void reCalculatePdMicromodulus(PD_Particles &particles, int dim) {
  arma::mat &data = particles.data();
  const ivec &idToCol = particles.getIdToCol_v();
  const ivec &colToId = particles.colToId();
  const int indexVolume = particles.getParamId("volume");
  const int indexMicromodulus = particles.getParamId("micromodulus");
  //   const int indexVolumeScaling = particles.getPdParamId("volumeScaling");
  const int indexDr0 = particles.getPdParamId("dr0");
  const double dimScaling = 2. * pow(dim, 2.);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int pId = colToId(i);
    double dRvolume = 0;

    const vector<pair<int, vector<double>>> &PDconnections =
        particles.pdConnections(pId);
    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int col_j = idToCol[id_j];
      //            double volumeScaling = con.second[indexVolumeScaling];
      const double dr0Len = con.second[indexDr0];

      dRvolume += dr0Len * data(col_j, indexVolume);
      //            dRvolume += dr0Len*data(col_j, indexVolume)*volumeScaling;
    }
    data(i, indexMicromodulus) *= dimScaling / dRvolume;
  }
}
//------------------------------------------------------------------------------
void reCalculatePdFractureCriterion(PD_Particles &particles, double G0,
                                    double delta, double h) {
  arma::mat &data = particles.data();
  const ivec &colToId = particles.colToId();
  const int indexMicromodulus = particles.getParamId("micromodulus");
  int indexS0 = particles.getParamId("s0");
  const double delta4 = delta * delta * delta * delta;
  const double delta5 = delta4 * delta;

  const ivec &idToCol = particles.getIdToCol_v();
  int indexVolume = particles.getParamId("volume");
  int indexVolumeScaling = particles.getPdParamId("volumeScaling");

  if (h != -1) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const double c = data(i, indexMicromodulus);
      const double s0 = sqrt(4 * G0 / (c * h * delta4));
      data(i, indexS0) = s0;
    }
  } else {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const int pId = colToId(i);

      const double c = data(i, indexMicromodulus);
      const double s0 = sqrt(10. * G0 / (c * M_PI * delta5));
      ;
      data(i, indexS0) = s0;

      double v = 0;

      vector<pair<int, vector<double>>> &PDconnections =
          particles.pdConnections(pId);
      for (auto &con : PDconnections) {
        int id_j = con.first;
        int col_j = idToCol[id_j];
        double volumeScaling = con.second[indexVolumeScaling];
        v += data(col_j, indexVolume) * volumeScaling;
      }
      double l_delta = pow((3. * v / (4. * M_PI)), 1. / 3.);
      double d5 = pow(l_delta, 5);
      double s_i = sqrt(10. * G0 / (c * M_PI * d5));
      //            double s_i = sqrt(60.*G0/(c*M_PI*delta5));
      data(i, indexS0) = s_i;
    }
  }
}
//------------------------------------------------------------------------------
void calculateRadius(PD_Particles &particles, int dim, double h) {
  arma::mat &data = particles.data();
  const int indexVolume = particles.getParamId("volume");
  const int indexRadius = particles.getParamId("radius");

  if (dim == 3) {
    const double scale = 0.9; // Shoudl be 0.704? Optimal sphere packing
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const double V = data(i, indexVolume);
      const double r = pow(3. * V / (4. * M_PI), 1. / 3.);
      data(i, indexRadius) = scale * r;
    }
  } else if (dim == 2) {
    const double scale = 0.9069;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const double V = data(i, indexVolume);
      const double r = sqrt(V / (M_PI * h));
      data(i, indexRadius) = scale * r;
      //            cout << data(i, indexRadius) << endl;
    }
  } else if (dim == 1) {
    const double scale = 1.0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const double V = data(i, indexVolume);
      const double r = 0.5 * V / (h * h);
      data(i, indexRadius) = scale * r;
    }
  } else {
    cerr << "Dimension not supported in radius calulations: " << dim << endl;
    throw;
  }
}
//------------------------------------------------------------------------------
void surfaceCorrection(PD_Particles &particles, vector<Force *> &forces,
                       double E, double nu, int dim) {
  const double strain = 0.001;
  vec3 scaleFactor;
  arma::mat &r = particles.r();
  const ivec &idToCol = particles.getIdToCol_v();
  const ivec &colToId = particles.colToId();
  arma::mat g = zeros(particles.nParticles(), dim);

  const int indexDr0 = particles.getPdParamId("dr0");
  const int indexForceScaling = particles.getPdParamId("forceScalingBond");
  const unsigned int nParticlesAndGhosts =
      particles.nParticles() + particles.nGhostParticles();
  // Stretching all particle in the x-direction
  scaleFactor(0) = strain;
  scaleFactor(1) = -nu * strain;
  scaleFactor(2) = -nu * strain;
  double W_infty = 0;

  for (int a = 0; a < dim; a++) {
    if (a == 1)
      scaleFactor.swap_rows(0, 1);
    else if (a == 2)
      scaleFactor.swap_rows(1, 2);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Loading the geometry
    for (unsigned int i = 0; i < nParticlesAndGhosts; i++) {
      const int col_i = i;

      for (int d = 0; d < dim; d++) {
        r(col_i, d) = (1 + scaleFactor(d)) * r(col_i, d);
      }
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Calculating the elastic energy density
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      const int id = colToId(i);
      double W = 0;

      for (Force *force : forces) {
        W += force->calculatePotentialEnergyDensity(id, i);
      }
      g(i, a) = W;
    }
    double median_g = arma::median(g.row(a));
    W_infty += median_g;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // Resetting the positions
    for (unsigned int i = 0; i < nParticlesAndGhosts; i++) {
      int col_i = i;

      for (int d = 0; d < dim; d++) {
        r(col_i, d) = r(col_i, d) / (1 + scaleFactor(d));
      }
    }
  }
  W_infty /= (dim);
  W_infty = 0.5 * E * pow(strain, 2);

  // Scaling the energy with the median energy, which we assume
  // to be the bulk energy
  for (int a = 0; a < dim; a++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < particles.nParticles(); i++) {
      int col_i = i;
      double g_i = g(col_i, a);
      double W = W_infty / g_i;
      g(col_i, a) = W;
    }
  }

// Calculating the scaling
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int pId = colToId(i);
    const int col_i = i;

    vector<pair<int, vector<double>>> &PDconnections =
        particles.pdConnections(pId);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int col_j = idToCol[id_j];

      const double dr0Len = con.second[indexDr0];
      vec3 n = (r.row(col_i).t() - r.row(col_j).t()) / dr0Len;

      vec3 g_mean;
      double G = 0;
      for (int d = 0; d < dim; d++) {
        g_mean(d) = 0.5 * (g(col_i, d) + g(col_j, d));
        G += pow(n(d) / g_mean(d), 2);
      }

      G = pow(G, -0.5);
      con.second[indexForceScaling] *= G;
    }
  }
}
//------------------------------------------------------------------------------
void applyInitialStrainStrain(PD_Particles &particles, double strain, int axis,
                              pair<double, double> area) {
  arma::mat &r = particles.r();
  const arma::mat &r0 = particles.r0();
  const unsigned int nParticlesAndGhosts =
      particles.nParticles() + particles.nGhostParticles();
  int counter = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < nParticlesAndGhosts; i++) {
    const int col_i = i;
    if (area.first <= r0(col_i, axis) && r0(col_i, axis) <= area.second) {
      const double rNew = (1 + strain) * r0(col_i, axis);
      r(col_i, axis) = rNew;
      counter++;
    }
  }
}
//------------------------------------------------------------------------------
void setPD_N3L(PD_Particles &particles) {
#ifdef USE_N3L
  const int indexMyPosistion = particles.registerPdParameter("myPosistion", -1);
  const ivec &colToId = particles.colToId();
  const ivec &idToCol = particles.getIdToCol_v();
  int n = 0;
  const arma::mat &r = particles.r();
  const int nParticles = particles.nParticles();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const int id_i = colToId(i);
    vector<pair<int, vector<double>>> &PDconnections_i =
        particles.pdConnections(id_i);

    for (auto &con_i : PDconnections_i) {
      const int id_j = con_i.first;
      const int j = idToCol[id_j];

      const vector<pair<int, vector<double>>> &PDconnections_j =
          particles.pdConnections(id_j);

      int found = false;
      int counter = 0;
      for (auto &con_j : PDconnections_j) {
        if (con_j.first == id_i) {
          found = true;
          break;
        }
        counter++;
      }

      con_i.second[indexMyPosistion] = counter;

      if (!found) {
        double dx, dy, dz;
        cerr << "N3L connections settings not found for" << id_i << " and "
             << id_j << endl;
        const vec &r_i = r.row(i).t();
        const vec &r_j = r.row(j).t();
        dx = r_i(0) - r_j(0);
        dy = r_i(1) - r_j(1);
        dz = r_i(2) - r_j(2);
        double dr = sqrt(dx * dx + dy * dy + dz * dz);
        cerr << dr << endl;
        n++;
      }
    }
  }
#endif
}
//------------------------------------------------------------------------------
void removeVoidConnections(PD_Particles &particles, Grid &grid,
                           const double delta, const double lc) {
  (void)grid;
  (void)delta;
  (void)lc;

  arma::mat &data = particles.data();
  ivec &idToCol = particles.getIdToCol_v();
  const ivec &colToId = particles.colToId();

  const mat &R0 = particles.r0();
  const int indexDr0 = particles.getPdParamId("dr0");
  const int indexRadius = particles.getParamId("radius");
  const int dim = grid.dim();
  //    const int m = 7; // Number of sample points
  const int m = 9; // Number of sample points
  const int indexConnected = particles.getPdParamId("connected");
  double sf = 1.0;

  if (dim == 3)
    sf = 1.35;
  //    sf = 1.8;
  if (dim == 2)
    sf = 1.55;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int pId_i = colToId.at(i);
    vector<pair<int, vector<double>>> &PDconnections =
        particles.pdConnections(pId_i);
    const vec &r_i = R0.row(i).t();
    ivec filled(m);
    arma::mat filledCenters = arma::zeros(M_DIM, m);
    const double radius_i = sf * data(i, indexRadius);

    for (auto &con : PDconnections) {
      const int id_j = con.first;
      const int col_j = idToCol[id_j];
      const double dr0_ij = con.second[indexDr0];
      const vec &r_j = R0.row(col_j).t();

      vec start = r_i;
      vec n_ij = (r_j - r_i) / dr0_ij;

      const double radius_j = sf * data(col_j, indexRadius);
      const double radius_ij = radius_i + radius_j;

      double r0 = radius_i;
      if (r_i(0) > r_j(0)) {
        start = r_j;
        n_ij = -n_ij;
        r0 = radius_j;
      }

      if (radius_ij >= dr0_ij)
        continue;

      const double spacing = (dr0_ij - radius_ij) / m;
      const double radius_spacing = 0.5 * spacing;

      for (int k = 0; k < m; k++) {
        filled(k) = 0;
        filledCenters.col(k) =
            start + (r0 + radius_spacing + k * spacing) * n_ij;
      }

      int nFilled = 0;

      for (auto &con_b : PDconnections) {
        const int id_b = con_b.first;
        const int b = idToCol[id_b];
        if (col_j == b)
          continue;
        const vec &r_b = R0.row(b).t();
        const double radius = sf * data(b, indexRadius);

        // Distance to each point
        for (int k = 0; k < m; k++) {
          const vec &checkPoint = filledCenters.col(k);
          const arma::vec &diff = checkPoint - r_b;
          double len_sq = 0;
          for (int d = 0; d < dim; d++) {
            len_sq += pow(diff(d), 2);
          }

          if (len_sq <= pow(radius + radius_spacing, 2)) {
            if (filled(k) < 1) {
              nFilled++;
            }
            filled(k) = 1;
          }
        }

        if (nFilled >= m)
          continue;
      }

      bool remove = false;

      for (int l = 0; l < m; l++) {
        if (filled(l) == 0) {
          remove = true;
        }
      }

      if (remove) {
        con.second[indexConnected] = 0;
      }
    }
  }

  // Enforcing symmetry
  int nFound = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int id_i = colToId.at(i);
    const vector<pair<int, vector<double>>> &PDconnections_i =
        particles.pdConnections(id_i);

    for (auto &con_j : PDconnections_i) {
      const int id_j = con_j.first;
      const bool connected = con_j.second[indexConnected];
      if (!connected) {
        vector<pair<int, vector<double>>> &PDconnections_j =
            particles.pdConnections(id_j);
        for (auto &con_k : PDconnections_j) {
          if (con_k.first == id_i) {
            con_k.second[indexConnected] = 0;
            nFound++;
          }
        }
      }
    }
  }

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  cout << "Void bonds removed due to symmetry = " << nFound << endl;
}
//------------------------------------------------------------------------------
void cleanUpPdConnections(PD_Particles &particles) {
  const int minNumberOfConnections = 5;
  const int indexConnected = particles.getPdParamId("connected");
  const ivec &colToId = particles.colToId();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int id_i = colToId(i);
    vector<pair<int, vector<double>>> &PDconnections_i =
        particles.pdConnections(id_i);
    vector<pair<int, vector<double>>> removeConnections;

    for (auto &con_j : PDconnections_i) {
      const bool connected = con_j.second[indexConnected];
      if (!connected) {
        removeConnections.push_back(con_j);
      }
    }

    for (auto &con_j : removeConnections) {
      PDconnections_i.erase(
          std::remove(PDconnections_i.begin(), PDconnections_i.end(), con_j),
          PDconnections_i.end());
    }

    if (PDconnections_i.size() <= minNumberOfConnections) {
      std::cerr << "Warning: node_id " << id_i
                << "\t\t has PDconnections.size(): " << PDconnections_i.size()
                << endl;
    }
  }
}
//------------------------------------------------------------------------------
void addFractures(PD_Particles &particles,
                  const vector<pair<double, double>> &domain,
                  const vector<double> &fracture) {
  //    double error_thres = 2.2204e-016;
  double error_thres = 1.e-10;
  cout << "Removing explicit fractures" << endl;
  const ivec &idToCol = particles.getIdToCol_v();
  const ivec &colToId = particles.colToId();
  const mat &R0 = particles.r0();
  const int indexConnected = particles.getPdParamId("connected");
  vec r2 = {fracture[0], fracture[2]};
  vec r3 = {fracture[1], fracture[3]};

  int broken_bonds = 0;
  for (unsigned int i = 0; i < particles.nParticles(); i++) {
    const int id_i = colToId(i);
    vector<pair<int, vector<double>>> &PDconnections =
        particles.pdConnections(id_i);
    vec r0 = R0.row(i).t();

    for (auto &con : PDconnections) {
      bool inverted = false;
      double x, y;

      const int id_j = con.first;
      const int j = idToCol[id_j];
      vec r1 = R0.row(j).t();

      // Special case where the particles are on top of each other in 3d
      if (r0[0] == r1[0] && r0[1] == r1[1])
        continue;

      if (abs(r3[0] - r2[0]) < error_thres &&
          abs(r1[1] - r0[1]) < error_thres) {
        x = r3[0];
        y = r1[1];
      } else if (abs(r3[1] - r2[1]) < error_thres &&
                 abs(r1[0] - r0[0]) < error_thres) {
        x = r1[0];
        y = r3[1];
      } else {
        // Special case where the x-axis is parallel with the axis cross
        if (r1[0] - r0[0] == 0 or r3[0] - r2[0] == 0) {
          inverted = true;
          r0 = {r0[1], r0[0]};
          r1 = {r1[1], r1[0]};
          r2 = {r2[1], r2[0]};
          r3 = {r3[1], r3[0]};
        }
        double a = (r1[1] - r0[1]) / (r1[0] - r0[0]);
        double b = r0[1] - r0[0] * a;

        double am = (r3[1] - r2[1]) / (r3[0] - r2[0]);
        double bm = r2[1] - r2[0] * am;

        if (a == am)
          continue;

        x = (bm - b) / (a - am);
        y = a * x + b;
      }

      vector<double> r = {x, y};

      bool intersect = true;
      for (int d = 0; d < 2; d++) {
        if (r[d] < std::min(r0[d], r1[d])) {
          if (std::abs(r[d] - std::min(r0[d], r1[d])) > error_thres)
            intersect = false;
        }
        if (r[d] > std::max(r0[d], r1[d])) {
          if (std::abs(r[d] - std::max(r0[d], r1[d])) > error_thres)
            intersect = false;
        }
        if (r[d] < std::min(r2[d], r3[d])) {
          if (std::abs(r[d] - std::min(r2[d], r3[d])) > error_thres)
            intersect = false;
        }
        if (r[d] > std::max(r2[d], r3[d])) {
          if (std::abs(r[d] - std::max(r2[d], r3[d])) > error_thres)
            intersect = false;
        }
      }
      if (inverted) {
        r0 = {r0[1], r0[0]};
        r1 = {r1[1], r1[0]};
        r2 = {r2[1], r2[0]};
        r3 = {r3[1], r3[0]};
      }

      if (intersect) {
        con.second[indexConnected] = 0;
        broken_bonds++;
      }
    }
  }

  cout << "broken_bonds: " << broken_bonds << endl;
}

//------------------------------------------------------------------------------
vector<array<int, 3>> possibleConfigurations2d(const int n) {
  vector<array<int, 3>> configurations;

  for (int nx = 1; nx <= n; nx++) {
    if (n % nx > 0)
      continue;
    const int ny = n / nx;
    configurations.push_back(array<int, 3>{nx, ny, 1});
  }

  return configurations;
}
//------------------------------------------------------------------------------
vector<array<int, 3>> possibleConfigurations3d(const int n) {
  vector<array<int, 3>> configurations;

  for (int nx = 1; nx <= n; nx++) {
    if (n % nx > 0)
      continue;

    const int nyz = n / nx;
    for (int ny = 1; ny <= n; ny++) {
      if (nyz % ny > 0)
        continue;

      const int nz = nyz / ny;
      configurations.push_back(array<int, 3>{nx, ny, nz});
    }
  }

  return configurations;
}
//------------------------------------------------------------------------------
vector<int> optimalConfigurationCores(const int nCores,
                                      const vector<double> &domain,
                                      const int dim) {

  vector<int> optimal = {1, 1, 1}; // default
  vector<double> nOptimal = {0, 0, 0};
  double totalLength = 0.;
  for (int d = 0; d < dim; d++) {
    nOptimal[d] = domain[d];
    totalLength += domain[d];
  }
  for (int d = 0; d < dim; d++) {
    nOptimal[d] = nCores * nOptimal[d] / totalLength;
  }

  vector<std::array<int, 3>> coreConfigurations;

  if (dim == 2)
    coreConfigurations = possibleConfigurations2d(nCores);
  else if (dim == 3)
    coreConfigurations = possibleConfigurations3d(nCores);

  double lowSigma = std::numeric_limits<double>::max();
  int optimalConfiguration = 0;
  int counter = 0;

  for (auto n : coreConfigurations) {
    double sigma = 0;
    for (int d = 0; d < dim; d++) {
      sigma += pow(n[d] - nOptimal[d], 2);
    }
    sigma = sqrt(sigma);
    if (sigma < lowSigma) {
      lowSigma = sigma;
      optimalConfiguration = counter;
    }
    counter++;
  }

  for (int d = 0; d < dim; d++) {
    optimal[d] = coreConfigurations[optimalConfiguration][d];
  }

#if 0
    for(int d=0; d<dim; d++)
    {
        cout << nOptimal[d] << " ";
    }
    cout << endl;

    for(auto n:coreConfigurations)
    {
        double sigma = 0;
        for(int d=0; d<dim; d++)
        {
            sigma += pow(n[d] - nOptimal[d], 2);
            cout << n[d] << "\t";
        }
        sigma = sqrt(sigma);
        cout << " \t s:" << sigma << endl;
    }
#endif
  return optimal;
}
//------------------------------------------------------------------------------
string getFileEnding(string filename) {
  std::vector<std::string> elements;
  boost::split(elements, filename, boost::is_any_of("."),
               boost::token_compress_on);
  return elements[elements.size() - 1];
}
//------------------------------------------------------------------------------
}
