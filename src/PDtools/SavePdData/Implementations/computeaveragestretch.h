#ifndef COMPUTEAVERAGESTRETCH_H
#define COMPUTEAVERAGESTRETCH_H

#include "PDtools/SavePdData/savepddata.h"
namespace PDtools
{

//------------------------------------------------------------------------------
class ComputeAverageStretch: public ComputeProperty
{
public:
    ComputeAverageStretch(PD_Particles &particles);
    ~ComputeAverageStretch();

    virtual void update(const pair<int, int> &pIdcol);
    virtual void init(const pair<int, int> &pIdcol);
private:
    int m_indexStretch;
    int m_indexAverageStretch;
    arma::mat *m_data;
};
//------------------------------------------------------------------------------
}

#endif // COMPUTEAVERAGESTRETCH_H
