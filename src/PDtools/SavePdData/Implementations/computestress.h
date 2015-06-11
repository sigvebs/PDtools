#ifndef COMPUTESTRESS_H
#define COMPUTESTRESS_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools
{
class Force;

//------------------------------------------------------------------------------
class ComputeStress: public ComputeProperty
{
public:
    ComputeStress(PD_Particles &particles, vector<Force*> &forces);
    ~ComputeStress();
    virtual void update(const pair<int, int> &pIdcol);
private:
    vector<Force *> &m_forces;
    arma::mat &m_data;
    int m_indexStress[6];
};
//------------------------------------------------------------------------------
}

#endif // COMPUTESTRESS_H
