#include <iostream>
#include <armadillo>
#include <memory>

#include "pdsolver.h"

int main(int argc, char** argv)
{
    using namespace arma;
    using namespace std;

    wall_clock timer;
    int myRank, nNodes;
    myRank = 0;
    nNodes = 1;

    if(argc != 2)
    {
        cerr << "usage: peridyn [path to config file]" << endl;
        return 1;
    }

    string configPath = argv[1];
    cout << configPath << endl;

    PdSolver Solver(configPath);
    Solver.initialize();

    timer.tic();
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

     return 0;
}
