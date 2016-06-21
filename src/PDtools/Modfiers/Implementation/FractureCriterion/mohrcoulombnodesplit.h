#ifndef MOHRCOULOMBNODESPLIT_H
#define MOHRCOULOMBNODESPLIT_H


#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulombNodeSplit : public Modifier
{
public:
    MohrCoulombNodeSplit(double mu, double C, double T);

    virtual void registerParticleParameters();
    virtual void initialize();
    virtual void evaluateStepOne();
    virtual void evaluateStepOne(const int id_i, const int i);
    virtual void evaluateStepTwo(const int id_i, const int i);
    virtual void evaluateStepTwo();

private:
    double m_C;
    double m_T;
    double m_d;
    double m_phi;
    double m_weight1;
    double m_weight2;
    arma::mat *m_data;
    std::unordered_map<int, int> * m_idToCol;

    int m_indexStress[6];
    int m_indexUnbreakable;
    int m_indexConnected;
    int m_indexCompute;
    int m_indexNewConnectionId_1;
    int m_indexNewConnectionId_2;
    int m_indexNormal[3];
    int m_indexBroken;
    int m_indexRadius;
    int m_indexVolume;
    int m_indexDr0;
    int m_indexBrokenNow;
    bool m_broken;
    vector<int> m_toBeDeleted;

    void copyParticleTo(int id_to, int old_col, int new_col);
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBNODESPLIT_H
