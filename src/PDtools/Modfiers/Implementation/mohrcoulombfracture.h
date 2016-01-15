#ifndef MOHRCOULOMBFRACTURE_H
#define MOHRCOULOMBFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulombFracture : public Modifier
{
public:
    MohrCoulombFracture(double mu, double C, double T, int dim);
    ~MohrCoulombFracture();

    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);

private:
    double m_C;
    double m_T;
    double m_d;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_dim;
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBFRACTURE_H
