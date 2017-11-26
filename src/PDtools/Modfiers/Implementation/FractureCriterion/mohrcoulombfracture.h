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
    MohrCoulombFracture(double mu, double C, double T);
    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);

private:
    double m_S0;
    double m_T;
    double m_d;
    double m_phi;
    double m_cos_theta;
    double m_sin_theta;
    double m_tan_theta;
    double m_tan2_theta;
    double m_k;
    double m_C0;

    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexBrokenNow;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBFRACTURE_H
