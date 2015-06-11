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
    ~PmbFracture();

    virtual void initialize();
    virtual void evaluateStepOne(const pair<int, int> &id_col);
    virtual void evaluateStepTwo(const pair<int, int> &id_col);

private:
    double m_alpha;
    int m_indexS0;
    int m_indexStretch;
    int m_indexUnbreakable;
    int m_indexS00;
    int m_indexS_tmp;
    arma::mat * m_data;
    std::unordered_map<int, int> * m_pIds;
};
//------------------------------------------------------------------------------
}

#endif // PMBFRACTURE_H
