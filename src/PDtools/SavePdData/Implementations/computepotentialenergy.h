#ifndef COMPUTEPOTENTIALENERGY_H
#define COMPUTEPOTENTIALENERGY_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools
{
class Force;

//------------------------------------------------------------------------------
class ComputePotentialEnergy: public ComputeProperty
{
public:
    ComputePotentialEnergy(PD_Particles &particles, vector<Force*> &forces);
    ~ComputePotentialEnergy();

    virtual void update(const pair<int, int> &pIdcol);
private:
    vector<Force *> &m_forces;
    arma::mat &m_data;
    int m_indexPotential;
};
//------------------------------------------------------------------------------
}

#endif // COMPUTEPOTENTIALENERGY_H
