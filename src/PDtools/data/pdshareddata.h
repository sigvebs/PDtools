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
//------------------------------------------------------------------------------
class PDsharedData
{
protected:
    map<string, double> m_data;
    int m_iteration;
    int m_nIterations;
    double m_t;
    double m_timeStep;
    pair<double, double> m_timeInterval;

    Domain *m_domain;
    Grid *m_pdGrid;
    Particles *m_particles;
    // or myParticles

    // MPI infomration
    int m_cId;
    int m_nProcs;

public:

    PDsharedData();
    ~PDsharedData();

    int iteration() const
    {
        return m_iteration;
    }

    void iteration(int _iteration)
    {
        m_iteration = _iteration;
    }

    int nIterations() const
    {
        return m_nIterations;
    }

    void nIterations(int _nIterations)
    {
        m_nIterations = _nIterations;
    }

    double t() const
    {
        return m_t;
    }

    void t(double _t)
    {
        m_t = _t;
    }

    double timeStep() const
    {
        return m_timeStep;
    }

    void timeStep(double _timeStep)
    {
        m_timeStep = _timeStep;
    }

    const pair<double, double> & timeInterval() const
    {
        return m_timeInterval;
    }

    void timeInterval(pair<double, double> _timeInterval)
    {
        m_timeInterval = _timeInterval;
    }

    Domain & domain() const
    {
        return *m_domain;
    }

    void domain(Domain *_domain)
    {
        m_domain= _domain;
    }

    Grid & pdGrid() const
    {
        return *m_pdGrid;
    }

    void pdGrid(Grid *_pdGrid)
    {
        m_pdGrid= _pdGrid;
    }

    Particles & particles() const
    {
        return *m_particles;
    }

    void particles(Particles *_particles)
    {
        m_particles= _particles;
    }

    int cID() const
    {
        return m_cId;
    }

    void cId(int _cId)
    {
        m_cId = _cId;
    }

    void registerNew(const string &key, double value);
    void updateData(const string &key, double value);
};
//------------------------------------------------------------------------------
}
#endif // PDSHAREDDATA_H
