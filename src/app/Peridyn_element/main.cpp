#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <armadillo>
#include <memory>
#include "pdsolver_e.h"

int main(int argc, char** argv)
{
    using namespace arma;
    using namespace std;
    wall_clock timer;
#ifdef USE_MPI
    MPI::Init (argc, argv);
    const int myRank = MPI::COMM_WORLD.Get_rank( );
    const int nCores = MPI::COMM_WORLD.Get_size( );
#else
    const int myRank = 0;
    const int nCores = 1;
#endif
    string configPath = argv[1];

    if(myRank == 0)
    {
        cout << "Starting with " << nCores << endl;
        cout << "input file: " << configPath << endl;
    }
    PdSolver_e Solver(configPath, myRank, nCores);
    timer.tic();
    Solver.initialize();

    Solver.solve();

    if(myRank == 0){
        double nSec = timer.toc();

        int hours = nSec/3600;
        int minutes = (nSec - hours*3600)/60;
        int seconds = (nSec - hours*3600 - minutes*60);

        cout << "Computed in "
             << hours << "h "
             << minutes << "m "
             << seconds << "s "
             << ", \t seconds: " << nSec
             << endl;
    }
#ifdef USE_MPI
    MPI::Finalize( );
#endif
    return 0;
}
