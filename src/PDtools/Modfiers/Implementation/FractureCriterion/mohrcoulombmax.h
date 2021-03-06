#ifndef MOHRCOULOMBAVERAGE_H
#define MOHRCOULOMBAVERAGE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class MohrCoulombMax : public Modifier {
public:
  MohrCoulombMax(double mu, double C, double T);

  virtual void registerParticleParameters();
  virtual void initialize();
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
  ivec *m_idToCol;

  int m_indexStress[6];
  int m_indexUnbreakable;
  int m_indexConnected;
  int m_indexCompute;
  //    int m_indexStressCenter;
  int m_indexBroken;
  bool m_broken;
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBAVERAGE_H
