#ifndef PDSHAREDDATA_H
#define PDSHAREDDATA_H

#include <map>
#include <assert.h>
#include <memory>

#include "Domain/domain.h"
#include "Particles/particles.h"
#include "Grid/grid.h"

using namespace std;
//------------------------------------------------------------------------------
namespace PDtools
{
class PDsharedData
{
protected:
    map<string, double> m_data;

public:
    int iteration;
    int nIterations;
    double t;
    double timeStep;
    pair<double, double> timeInterval;

    Domain *domain;
    Grid *pdGrid;
    Particles *particles;
    // or myParticles

    // MPI infomration
    int cId;
    int nProcs;

    PDsharedData();
    ~PDsharedData();

    void registerNew(const string &key, double value);
    void updateData(const string &key, double value);};
//------------------------------------------------------------------------------
}
#endif // PDSHAREDDATA_H
