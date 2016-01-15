#ifndef MOHRCOULOMBBONDFRACTURE_H
#define MOHRCOULOMBBONDFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulombBondFracture : public Modifier
{
public:
    MohrCoulombBondFracture(double mu, double C, double T, int dim);
    ~MohrCoulombBondFracture();

    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);

protected:
    double m_C;
    double m_T;
    double m_d;
    arma::mat *m_data;
    arma::mat *m_r;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_dim;
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBBONDFRACTURE_H
