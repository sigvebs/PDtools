#ifndef PMBFRACTURE_H
#define PMBFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class PmbFracture: public Modifier
{
public:
    PmbFracture(double alpha);

    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);
    virtual void updateStepOne(const int id, const int i);

private:
    double m_alpha;
    int m_indexS0;
    int m_indexStretch;
    int m_indexUnbreakable;
    int m_indexS00;
    int m_indexConnected;
    int m_indexS_tmp;
    int m_indexBrokenNow;
    double m_s00;
    arma::mat * m_data;
    std::unordered_map<int, int> * m_idToCol;
};
//------------------------------------------------------------------------------
}

#endif // PMBFRACTURE_H
