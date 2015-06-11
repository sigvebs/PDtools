#ifndef ADRMOHRCOULOMBFRACTURE_H
#define ADRMOHRCOULOMBFRACTURE_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
class Force;
//------------------------------------------------------------------------------
class ADRmohrCoulombFracture : public Modifier
{
public:
    ADRmohrCoulombFracture(double mu, double C, double T);
    ~ADRmohrCoulombFracture();

    virtual void initialize();
    virtual void evaluateStepOne(const pair<int, int> &pIdcol);
    virtual void evaluateStepTwo(const pair<int, int> &pIdcol);
    virtual void evaluateStepTwo();

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
    pair<int, pair<int, vector<double>> *> m_maxPId;
    double m_maxStress;
};
//------------------------------------------------------------------------------
}
#endif // ADRMOHRCOULOMBFRACTURE_H
