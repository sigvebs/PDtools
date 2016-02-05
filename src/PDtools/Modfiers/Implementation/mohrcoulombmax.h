#ifndef MOHRCOULOMBAVERAGE_H
#define MOHRCOULOMBAVERAGE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulombMax : public Modifier
{
public:
    MohrCoulombMax(double mu, double C, double T, int dim);

    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);
    virtual void evaluateStepTwo(const int id_i, const int i);
    virtual void evaluateStepTwo();

private:
    double m_C;
    double m_T;
    double m_d;
    double m_phi;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_dim;
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
    int m_indexStressCenter;
    int m_indexBroken;
    bool m_broken;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBAVERAGE_H
