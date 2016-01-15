#ifndef MICROPOLARFRACTURE_H
#define MICROPOLARFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------

class MicropolarFracture : public Modifier
{
public:
    MicropolarFracture(double theta);
    ~MicropolarFracture();

    virtual void
    registerParticleParameters();
    virtual void
    initialize();
    virtual void
    evaluateStepOne(const int id_i, const int i);
private:
    double m_thetaCritical;
    double m_s00;
    int m_indexS0;
    int m_indexStretch;
    int m_indexUnbreakable;
    int m_indexS00;
    int m_indexConnected;
    int m_indexTheta;
    arma::mat * m_data;
    std::unordered_map<int, int> * m_idToCol;
};
//------------------------------------------------------------------------------
}
#endif // MICROPOLARFRACTURE_H
