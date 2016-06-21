#ifndef PDSOLVER_H
#define PDSOLVER_H

#include <PDtools.h>
#include <libconfig.h++>

using namespace std;

//------------------------------------------------------------------------------
class PdSolver_e
{
public:
    PdSolver_e(string configPath, int myRank, int nMpiNodes);
    ~PdSolver_e();

    int initialize();
    void solve();

protected:
    const string m_configPath;
    const int m_myRank;
    const int m_nCores;
    bool isRoot = 0;
    libconfig::Config m_cfg;
    PDtools::PD_Particles m_discretization;
    PDtools::Grid m_grid;
    PDtools::Solver *solver;
};
//------------------------------------------------------------------------------
#endif // PDSOLVER_H
