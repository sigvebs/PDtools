#include "adr.h"

#if USE_MPI
#include <mpi.h>
#endif

#include "PDtools/Force/force.h"
#include "PDtools/Grid/grid.h"
#include "PDtools/Modfiers/modifier.h"
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/PdFunctions/pdfunctions.h"
#include "PDtools/SavePdData/savepddata.h"

#include <limits> // TMP
using namespace arma;

namespace PDtools
//------------------------------------------------------------------------------
{
ADR::ADR() {
  m_dt = 1.0;
  m_c = 0;
  m_du_u = 0;
}
//------------------------------------------------------------------------------
void ADR::solve() {
  initialize();
  checkInitialization();
  save(0);

  // Looping over all time, particles and components.
  for (int i = 0; i < m_steps; i++) {
    stepForward(i);
  }
}
//------------------------------------------------------------------------------
void ADR::checkInitialization() {}
//------------------------------------------------------------------------------
void ADR::initialize() {
  mat &F = m_particles->F();
  mat &v = m_particles->v();
  mat &Fold = m_particles->Fold();
  mat &r_prev = m_particles->r_prev();
  const mat &r0 = m_particles->r0();

  const int nParticles = m_particles->nParticles();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    for (int d = 0; d < m_dim; d++) {
      F(i, d) = 0;
      Fold(i, d) = 0;
      v(i, d) = 0;
      r_prev(i, d) = r0(i, d);
    }
  }

  updateGridAndCommunication();
  //-----------------------------
  calculateForces(0);
  calculateStableMass();
  updateProperties(0);
  Solver::initialize();
}
//------------------------------------------------------------------------------
void ADR::stepForward(int i) {
  modifiersStepOne();
  //    updateGridAndCommunication(); // NEW - tmp
  applyBoundaryConditions();
  updateGridAndCommunication();

  m_globalError = 2 * m_errorThreshold;
  iterate(m_maxSteps);
  //    updateGridAndCommunication(); // NEW - tmp
  updateProperties(i);

  updateGridAndCommunication();
  modifiersStepTwo();
  updateGridAndCommunication();
  save(i + 1);
  m_t += m_dt;

  if (m_myRank == 0)
    cout << "i = " << i << " " << m_globalError << endl;
}
//------------------------------------------------------------------------------
void ADR::iterate(int maxNumberOfSteps) {
  //    int nFractureRelaxationSteps = m_maxStepsFracture;
  //    int nFractureRelaxationSteps = 2*minSteps;
  //    const ivec &colToId = m_particles->colToId();
  //    int minSteps = 500;
  //        int nRelaxSteps = 30;
  const int minSteps = m_minSteps;

  int counter = 0;
  int steps = 0;
  int relaxSteps = minSteps;
  //    maxNumberOfSteps = 1; // TMP
  // -------------------------------------------------------------------------

  for (counter = 0; counter < maxNumberOfSteps; counter++) {
    //        cout << counter << " global_error:" << m_globalError << "
    //        err_threshold:"<< m_errorThreshold << endl;
    integrateStepOne();
    updateGridAndCommunication();
    zeroForces();
    calculateForces(0);
    staticModifiers();
    integrateStepTwo();
    updateGridAndCommunication();

    //        if(m_myRank == 0)
    //            cout << "c:" << counter << " err:" << m_globalError <<  endl;
    /*
    if(relaxSteps == 0) {
        int continueState = 0;
        int nParticles = m_particles->nParticles();

        // It now loops over all properties, not only the ones needed in
        // static computations.
//            updateGridAndCommunication(); // NEW - tmp
        updateProperties(0);
        updateGhosts();
//            updateGridAndCommunication(); // NEW - tmp

        // Forces modifiers
        for(Force *oneBodyForce:m_oneBodyForces) {
            if(!oneBodyForce->getHasStaticModifier())
                continue;

            for(int i=0; i<nParticles; i++) {
                const int id = colToId(i);
                oneBodyForce->evaluateStatic(id, i);
            }
            if(oneBodyForce->getContinueState() > 0 )
                continueState = 1;
        }

//            updateGridAndCommunication(); // NEW - tmp

        for(Modifier *modifier:m_qsModifiers) {
            if(!modifier->hasStepOne())
                continue;

            for(int i=0; i<nParticles; i++) {
                const int id = colToId(i);
                modifier->evaluateStepOne(id, i);
            }
        }

        updateGhosts();
//            updateGridAndCommunication(); // NEW - tmp
        nParticles = m_particles->nParticles();

        for(Modifier *modifier:m_qsModifiers) {
            if(!modifier->hasStepTwo())
                continue;

            for(int i=0; i<nParticles; i++) {
                const int id = colToId(i);
                modifier->evaluateStepTwo(id, i);
            }
        }

        for(Modifier *modifier:m_qsModifiers) {
            modifier->evaluateStepTwo();
        }
//            updateGridAndCommunication(); // NEW - tmp

        for(Modifier *modifier:m_qsModifiers) {
            if(modifier->state())
                continueState = 1;
        }
//            updateGridAndCommunication(); // NEW - tmp

#if USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &continueState, 1, MPI_INT, MPI_SUM,
MPI_COMM_WORLD);
#endif
        if(continueState > 0) {
            m_globalError = 2*m_errorThreshold;
            counter -= m_maxStepsFracture;
            relaxSteps = nRelaxSteps;
            if(counter < 0)
                counter = 0;
        }
    }
    */

    //        updateGridAndCommunication(); // NEW - tmp
    relaxSteps--;
    if (relaxSteps < 0)
      relaxSteps = 0;
    steps++;
    if ((m_globalError < m_errorThreshold) && (counter > minSteps))
      break;
  }

  if (m_myRank == 0)
    cout << "counter:" << steps << " err:" << m_globalError << endl;

  m_particles->uppdateR_prev();
  // -------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void ADR::integrateStepOne() {
  //    ivec & colToId = m_particles->colToId();
  //    const arma::imat & isStatic = m_particles->isStatic();
  using namespace arma;

  const double alpha = (2. * m_dt) / (2. + m_c * m_dt);
  const double beta = (2. - m_c * m_dt) / (2. + m_c * m_dt);

  const vec &stableMass = m_particles->stableMass();
  mat &r = m_particles->r();
  mat &v = m_particles->v();
  mat &F = m_particles->F();
  mat &Fold = m_particles->Fold();
  //    mat & r0 = m_particles->r0();

  const int nParticles = m_particles->nParticles();

  double maxU = 0;
  double avgU = 0;
  int np = 0;

  // TMP for convergence test
  const mat &r_prev = m_particles->r_prev();

  // Pointers
  double *x = r.colptr(0);
  double *y = r.colptr(1);
  double *z = r.colptr(2);
  const double *x_p = r_prev.colptr(0);
  const double *y_p = r_prev.colptr(1);
  const double *z_p = r_prev.colptr(2);
  double *vx = v.colptr(0);
  double *vy = v.colptr(1);
  double *vz = v.colptr(2);

  const double *Fx = F.colptr(0);
  const double *Fy = F.colptr(1);
  const double *Fz = F.colptr(2);
  double *Foldx = Fold.colptr(0);
  double *Foldy = Fold.colptr(1);
  double *Foldz = Fold.colptr(2);
  const double *sm = stableMass.colptr(0);

  if (m_dim == 3) {
    for (int i = 0; i < nParticles; i++) {
      //            if(isStatic[i])
      //                continue;

      const double sm_i = sm[i];
      vx[i] = beta * vx[i] + alpha * Fx[i] / sm_i;
      vy[i] = beta * vy[i] + alpha * Fy[i] / sm_i;
      vz[i] = beta * vz[i] + alpha * Fz[i] / sm_i;
      x[i] += vx[i] * m_dt;
      y[i] += vy[i] * m_dt;
      z[i] += vz[i] * m_dt;

      Foldx[i] = Fx[i];
      Foldy[i] = Fy[i];
      Foldz[i] = Fz[i];

      const double l_u =
          pow(x[i] - x_p[i], 2) + pow(y[i] - y_p[i], 2) + pow(z[i] - z_p[i], 2);

      if (l_u > 1.e-22) {
        const double sqrt_lu = sqrt(l_u);
        maxU = std::max(sqrt_lu, maxU);
        avgU += sqrt_lu;
        np++;
      }
    }
  } else { // dim 2
    for (int i = 0; i < nParticles; i++) {
      //            if(isStatic[i])
      //                continue;

      const double sm_i = sm[i];
      vx[i] = beta * vx[i] + alpha * Fx[i] / sm_i;
      vy[i] = beta * vy[i] + alpha * Fy[i] / sm_i;
      x[i] += vx[i] * m_dt;
      y[i] += vy[i] * m_dt;
      Foldx[i] = Fx[i];
      Foldy[i] = Fy[i];

      const double l_u = pow(x[i] - x_p[i], 2) + pow(y[i] - y_p[i], 2);

      if (l_u > 1.e-22) {
        const double sqrt_lu = sqrt(l_u);
        maxU = std::max(sqrt_lu, maxU);
        avgU += sqrt_lu;
        np++;
      }
    }
  }
#if USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxU, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &avgU, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  //    if(m_myRank == 0)
  //        cout << maxU << "\t" <<  avgU << "\t" <<  maxU/avgU << "\t" << np <<
  //        "\t" << m_globalError << endl;
  const double nTot = m_particles->totParticles();
  const double npTot = np / nTot;

  if (m_counter % 10 == 0 && npTot > 0.8) {
    avgU /= np;
    m_globalError = fabs(maxU / avgU - m_du_u);
    m_du_u = maxU / avgU;

    if (std::isnan(m_globalError))
      m_globalError = 2 * m_errorThreshold;

    //        if(m_myRank == 0)
    //            cout << maxU << "\t" <<  avgU << "\t" <<  maxU/avgU << "\t" <<
    //            m_globalError << endl;
  }
  m_counter++;
}
//------------------------------------------------------------------------------
void ADR::integrateStepTwo() {
  const mat &v = m_particles->v();
  const mat &F = m_particles->F();
  const mat &Fold = m_particles->Fold();
  const vec &stableMass = m_particles->stableMass();
  const arma::imat &isStatic = m_particles->isStatic();

  // Pointers
  const double *vx = v.colptr(0);
  const double *vy = v.colptr(1);
  const double *vz = v.colptr(2);
  const double *Fx = F.colptr(0);
  const double *Fy = F.colptr(1);
  const double *Fz = F.colptr(2);
  const double *Foldx = Fold.colptr(0);
  const double *Foldy = Fold.colptr(1);
  const double *Foldz = Fold.colptr(2);
  const double *sm = stableMass.colptr(0);

  // Calculating the damping coefficient
  double numerator = 0;
  double denominator = 0;
  int const nParticles = m_particles->nParticles();

  if (m_dim == 3) {
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : numerator, denominator)
#endif
    for (int i = 0; i < nParticles; i++) {
      if (isStatic[i])
        continue;

      const double sm_i = sm[i] * m_dt;
      numerator -= vx[i] * (Fx[i] - Foldx[i]) / sm_i;
      numerator -= vy[i] * (Fy[i] - Foldy[i]) / sm_i;
      numerator -= vz[i] * (Fz[i] - Foldz[i]) / sm_i;
      denominator += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
  } else { // dim 2
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : numerator, denominator)
#endif
    for (int i = 0; i < nParticles; i++) {
      if (isStatic(i))
        continue;

      const double sm_i = sm[i] * m_dt;
      numerator -= vx[i] * (Fx[i] - Foldx[i]) / sm_i;
      numerator -= vy[i] * (Fy[i] - Foldy[i]) / sm_i;
      denominator += vx[i] * vx[i] + vy[i] * vy[i];
    }
  }
#if USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &numerator, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &denominator, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  m_c = 0;
  if (denominator > 0) {
    if (numerator / denominator > 0) {
      m_c = 2 * sqrt(numerator / denominator);
    }
  }

  if (m_c >= 2.0) {
    m_c = 1.99;
  }
}
//------------------------------------------------------------------------------
void ADR::staticModifiers() {
  for (Modifier *modifier : m_boundaryModifiers) {
    modifier->staticEvaluation();
  }
}
//------------------------------------------------------------------------------
void ADR::calculateStableMass() {
#if USE_MPI
  m_mainGrid->clearGhostParticles();
  exchangeInitialGhostParticles(*m_mainGrid, *m_particles);
#endif
  const ivec &colToId = m_particles->colToId();
  arma::vec &stableMass = m_particles->stableMass();
  const int nParticles = m_particles->nParticles();

  double average_stableMass = 0;
  double min_stableMass = std::numeric_limits<double>::max();
  double max_stableMass = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const int id = colToId(i);
    double stableMassMax = 0;
    for (Force *oneBodyForce : m_oneBodyForces) {
      const double sm = oneBodyForce->calculateStableMass(id, i, m_dt);
      stableMassMax = std::max(stableMassMax, sm);
    }
    //        cout << "sm:" << stableMassMax << endl;
    stableMass(i) = stableMassMax;
    average_stableMass += stableMassMax;

    if (stableMassMax < min_stableMass)
      min_stableMass = stableMassMax;
    if (stableMassMax > max_stableMass)
      max_stableMass = stableMassMax;
  }
#if USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &average_stableMass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_stableMass, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max_stableMass, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  average_stableMass = average_stableMass / m_nCores / nParticles;

  if (m_myRank == 0) {
    cout << "average stable mass:" << average_stableMass << endl;
    cout << "min stable mass:" << min_stableMass << endl;
    cout << "max stable mass:" << max_stableMass << endl;
  }
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

// TMP fix for zero stable mass nodes
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nParticles; i++) {
    const int id = colToId(i);
    double sm = stableMass(i);

    if (sm < average_stableMass / 10) {
      cerr << "Warning: stable mass is close to zero for node " << id << endl;
      cerr << "Setting it to average_stableMass/10:" << average_stableMass / 10
           << endl;
      stableMass(i) = average_stableMass / 10;
    }
  }
}
//------------------------------------------------------------------------------
void ADR::updateGridAndCommunication() {
  m_mainGrid->clearParticles();
  updateGrid(*m_mainGrid, *m_particles, true);

#if USE_MPI
  exchangeGhostParticles(*m_mainGrid, *m_particles);

  // Updating lists in modifiers
  int counter = 0;
  for (Modifier *modifier : m_boundaryModifiers) {
    updateModifierLists(*modifier, *m_particles, counter);
    counter++;
  }
#endif

  //    updateElementQuadrature(*m_particles);
  //    if(m_myRank == 2) {
  //        mat & r = m_particles->r();
  //        unordered_map<int, int>  & idToCol = m_particles->idToCol();
  //        const int id_i = 14;
  //        const int col = idToCol[id_i];
  ////        updateGridAndCommunication();
  //        cout << "now: " << r(col, 1) << " c:" << col <<  endl;
  //    }
}
//------------------------------------------------------------------------------
}
