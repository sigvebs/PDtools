#ifndef MOHRCOULOMBFRACTURE_H
#define MOHRCOULOMBFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
class Force;
//------------------------------------------------------------------------------
class MohrCoulombFracture : public Modifier
{
public:
    MohrCoulombFracture(double mu, double C, double T);
    ~MohrCoulombFracture();

    virtual void initialize();
    virtual void evaluateStepOne(const pair<int, int> &pIdcol);
    virtual void evaluateStepTwo(const pair<int, int> &pIdcol);

    void addForce(Force *force);

private:
    double m_C;
    double m_T;
    double m_d;
    vector<Force *> m_forces;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_pIds;

    int m_indexStress[6];
    int m_dim;
    int m_indexUnbreakable;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBFRACTURE_H
