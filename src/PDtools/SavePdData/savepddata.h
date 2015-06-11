#ifndef SAVEPDDATA_H
#define SAVEPDDATA_H

#include <vector>
#include <iostream>
#include <armadillo>
#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Particles/saveparticles.h"
//#include <boost/filesystem.hpp>

//------------------------------------------------------------------------------
namespace PDtools
{
using namespace std;

// Forward declerations
class Force;

//------------------------------------------------------------------------------
class ComputeProperty
{
public:
    ComputeProperty(PD_Particles &particles);

    virtual void update(const pair<int, int> &pIdcol) = 0;
    virtual void init(const pair<int, int> &pIdcol);
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

    void initialize();
    void evaluate(double t, int i);

    void setSavePath(string savePath);
    void setParticles(PD_Particles *particles);
    void addParameter(std::string param);
    void setForces(std::vector<Force *> &oneBodyForces);

private:
    PD_Particles *m_particles = nullptr;
    std::vector<std::string> m_saveParameters = {"id", "x", "y"};
    std::vector<Force *> *m_oneBodyForces;
    std::vector<ComputeProperty *> m_computeProperties;
    SaveParticles *saveParticles;
    string m_savePath;
};
//------------------------------------------------------------------------------
}
#endif // SAVEPDDATA_H
