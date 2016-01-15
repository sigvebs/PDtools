#ifndef BONDENERGYFRACTURE_H
#define BONDENERGYFRACTURE_H


#include <unordered_map>
#include "PDtools/Force/force.h"
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class BondEnergyFracture : public Modifier
{
public:
    BondEnergyFracture(double delta, double G0, int dim, vector<Force*> &forces, double h=1.);
    ~BondEnergyFracture();

    virtual void initialize();
    virtual void evaluateStepOne(const int id_i, const int i);
    virtual void evaluateStepTwo(); protected:
    vector<Force*> &m_forces;
    double m_wc; // Critical energy density
    int m_indexStretch;
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexDr0;
    int m_broken;
    int m_indexForceScaling;
    int m_indexMicromodulus;
    double  m_G;
    arma::mat * m_data;
    std::unordered_map<int, int> * m_idToCol;
};
//------------------------------------------------------------------------------
}
#endif // BONDENERGYFRACTURE_H
