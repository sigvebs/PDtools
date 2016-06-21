#ifndef STRAINFRACTURE_H
#define STRAINFRACTURE_H


#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
class StrainFracture : public Modifier
{
public:
    StrainFracture(double Eeq, double Evol, int dim);

    virtual void
    registerParticleParameters();

    virtual void
    evaluateStepOne(const int id_i, const int i);

protected:
    double m_Eeq;
    double m_Evol;
    int m_dim;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
    int m_indexBrokenNow;
    int m_indexStrain[6];
};
//------------------------------------------------------------------------------
}

#endif // STRAINFRACTURE_H
