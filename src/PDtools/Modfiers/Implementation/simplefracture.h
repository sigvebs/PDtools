#ifndef SIMPLEFRACTURE_H
#define SIMPLEFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class SimpleFracture : public Modifier
{
public:
    SimpleFracture(double alpha);
    ~SimpleFracture();

    virtual void initialize();
    virtual void evaluateStepOne(const pair<int, int> &pIdcol);
    virtual void evaluateStepTwo();

private:
    double m_alpha;
    int m_indexStretch;
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexDr0;
    int m_indexVolume;
    int m_indexForceScaling;
    int m_indexMicromodulus;
    bool m_broken;
    arma::mat * m_data;
    std::unordered_map<int, int> * m_pIds;
};
//------------------------------------------------------------------------------
}
#endif // SIMPLEFRACTURE_H
