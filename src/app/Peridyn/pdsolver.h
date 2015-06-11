#ifndef PDSOLVER_H
#define PDSOLVER_H

#include <PDtools.h>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

//------------------------------------------------------------------------------
class PdSolver
{
public:
    PdSolver(string configPath);
    ~PdSolver();

    void initialize();
    void solve();

protected:
    string m_configPath;
    Config m_cfg;
    PDtools::PD_Particles m_particles;
    PDtools::Grid m_grid;
    PDtools::Solver *solver;
};
//------------------------------------------------------------------------------
#endif // PDSOLVER_H
