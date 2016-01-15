#ifndef SAVEPDDATA_H
#define SAVEPDDATA_H

#include <vector>
#include <iostream>
#include <armadillo>
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Particles/saveparticles.h"
//#include <boost/filesystem.hpp>
#ifdef USE_MPI
#include <mpi.h>
#endif

//------------------------------------------------------------------------------
namespace PDtools
{
using namespace std;

// Forward declerations
class Force;
class Grid;

//------------------------------------------------------------------------------
class ComputeProperty
{
public:
    ComputeProperty(PD_Particles &particles);

    virtual void update(const int id_i, const int i) = 0;
    virtual void init(const int id_i, const int i);
protected:
    PD_Particles &m_particles;
};
//------------------------------------------------------------------------------
class SavePdData
{
public:
    SavePdData();
    SavePdData(std::vector<std::string> saveParameters);
    ~SavePdData();

    void
    initialize();
    void
    evaluate(double t, int i);

    void
    saveData(double t, int i);

    void
    setSavePath(string savePath);
    void
    setParticles(PD_Particles *particles);
    void
    addParameter(std::string param);
    void
    setForces(std::vector<Force *> &oneBodyForces);
    void
    setScaling(const double E0, const double L0, const double v0, const double t0, const double rho0);
    void
    setRankAndCores(int rank, int cores);
    void
    setGrid(Grid *grid);
    bool writeBinary() const;
    void setWriteBinary(bool writeBinary);

    int
    updateFrquency() const;

    void
    setUpdateFrquency(int updateFrquency);

    std::vector<pair<string, int> >
    neededProperties() const;

    void
    setDim(int dim);

private:
    int m_updateFrquency;
    PD_Particles *m_particles = nullptr;
    std::vector<std::string> m_saveParameters = {"id", "x", "y"};
    std::vector<std::pair<std::string, double>> m_saveparam_scale;
    std::vector<Force *> *m_oneBodyForces;
    std::vector<ComputeProperty *> m_computeProperties;
    SaveParticles *saveParticles;
    string m_savePath;
    Grid *m_mainGrid;
    std::vector<pair<string, int>> m_neededProperties; // name and update frequency
    bool m_writeBinary = false;

    double m_E0 = 1.;
    double m_L0 = 1.;
    double m_v0 = 1.;
    double m_t0 = 1.;
    double m_rho0 = 1;

    int m_myRank = 0;
    int m_nCores = 1;
    int m_dim;
};
//------------------------------------------------------------------------------
}
#endif // SAVEPDDATA_H
