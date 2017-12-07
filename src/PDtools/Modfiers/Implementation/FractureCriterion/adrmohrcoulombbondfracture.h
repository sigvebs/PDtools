#ifndef ADRMOHRCOULOMBBONDFRACTURE_H
#define ADRMOHRCOULOMBBONDFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class AdrMohrCoulombBondFracture : public Modifier {
public:
  AdrMohrCoulombBondFracture(double mu, double C, double T);
  virtual void initialize();
  virtual void evaluateStepTwo(const int id_i, const int i);
  virtual void evaluateStepTwo();

private:
  double m_C;
  double m_T;
  double m_d;
  arma::mat *m_data;
  std::unordered_map<int, int> *m_idToCol;

  int m_indexStress[6];
  int m_indexUnbreakable;
  int m_indexConnected;
  pair<int, int> m_maxPId;
  double m_maxStress;
  int m_indexBrokenNow;
  int m_nStressElements;
};
//------------------------------------------------------------------------------
}
#endif // ADRMOHRCOULOMBBONDFRACTURE_H
