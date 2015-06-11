#ifndef COMPUTEDAMAGE_H
#define COMPUTEDAMAGE_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools
{

//------------------------------------------------------------------------------
class ComputeDamage: public ComputeProperty
{
public:
    ComputeDamage(PD_Particles &particles);
    ~ComputeDamage();

    virtual void update(const pair<int, int> &pIdcol);
    virtual void init(const pair<int, int> &pIdcol);
private:
    int m_indexDamage;
    int m_indexMaxPdConnections;
    arma::mat *m_data;
};
//------------------------------------------------------------------------------
}
#endif // COMPUTEDAMAGE_H
