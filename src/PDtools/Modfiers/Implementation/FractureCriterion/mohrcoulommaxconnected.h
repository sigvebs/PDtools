#ifndef MOHRCOULOMMAXCONNECTED_H
#define MOHRCOULOMMAXCONNECTED_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulomMaxConnected : public Modifier
{
public:
    MohrCoulomMaxConnected(double mu, double C, double T);
    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);
    virtual void evaluateStepTwo(const int id_i, const int i);

private:
    double m_C;
    double m_T;
    double m_d;
    double m_phi;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_indexNormal[3];
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
    int m_indexBrokenNow;
    bool m_broken;
    int m_indexBroken;
    int m_indexBrokenId;
    double m_cosTheta;
    double m_sinTheta;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMMAXCONNECTED_H
