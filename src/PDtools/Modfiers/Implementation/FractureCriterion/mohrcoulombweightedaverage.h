#ifndef MOHRCOULOMBWEIGHTEDAVERAGE_H
#define MOHRCOULOMBWEIGHTEDAVERAGE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class MohrCoulombWeightedAverage : public Modifier {
public:
  MohrCoulombWeightedAverage(double mu, double C, double T);

  virtual void registerParticleParameters();
  virtual void initialize();
  virtual void evaluateStepOne(const int id_i, const int i);

private:
  double m_C;
  double m_T;
  double m_d;
  double m_phi;
  arma::mat *m_data;
  std::unordered_map<int, int> *m_idToCol;
  double m_weights[2];

  int m_indexStress[6];
  int m_indexUnbreakable;
  int m_indexConnected;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBWEIGHTEDAVERAGE_H
