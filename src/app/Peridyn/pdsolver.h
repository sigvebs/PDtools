#ifndef PDSOLVER_H
#define PDSOLVER_H

#include <PDtools.h>
#include <libconfig.h++>
#include <memory>
#include <string>

//------------------------------------------------------------------------------
class PdSolver {
public:
  PdSolver(std::string configPath, int myRank, int nMpiNodes);
  ~PdSolver();

  int initialize();
  void solve();

private:
  void setDomain();

  const std::string m_configPath;
  const int m_myRank;
  const int m_nCores;
  bool isRoot = 0;
  libconfig::Config m_cfg;
  PDtools::PD_Particles m_particles;
  PDtools::Grid m_grid;
  PDtools::Solver *solver;
};
//------------------------------------------------------------------------------
#endif // PDSOLVER_H
