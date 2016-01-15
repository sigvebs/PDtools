#ifndef COMPUTEKINETICENERGY_H
#define COMPUTEKINETICENERGY_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools
{

//------------------------------------------------------------------------------
class ComputeKineticEnergy: public ComputeProperty
{
public:
    ComputeKineticEnergy(PD_Particles &particles);
    ~ComputeKineticEnergy();

    virtual void update(const int id_i, const int i);

private:
    const int m_dim = 3;
    int m_indexKE;
    int m_indexVolume;
    int m_indexRho;
    arma::mat *m_v;
    arma::mat *m_data;
};
//------------------------------------------------------------------------------
}

#endif // COMPUTEKINETICENERGY_H
