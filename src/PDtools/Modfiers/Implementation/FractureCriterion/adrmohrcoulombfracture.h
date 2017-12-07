#ifndef ADRMOHRCOULOMBFRACTURE_H
#define ADRMOHRCOULOMBFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class ADRmohrCoulombFracture : public Modifier {
public:
  ADRmohrCoulombFracture(double mu, double C, double T);
  virtual void initialize();
  virtual void evaluateStepTwo(const int id_i, const int i);
  virtual void evaluateStepTwo();

private:
  double m_C;
  double m_T;
  double m_d;
  double m_phi;
  double m_cos_theta;
  double m_sin_theta;
  double m_tan_theta;
  double m_tan2_theta;

  arma::mat *m_data;
  ivec *m_idToCol;

  int m_indexStress[6];
  int m_indexUnbreakable;
  int m_indexConnected;
  int m_indexBrokenNow;
  pair<int, int> m_maxPId;
  double m_maxStress;
  int m_nStressElements;
};
//------------------------------------------------------------------------------
}
#endif // ADRMOHRCOULOMBFRACTURE_H
