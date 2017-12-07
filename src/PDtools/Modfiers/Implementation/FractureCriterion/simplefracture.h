#ifndef SIMPLEFRACTURE_H
#define SIMPLEFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class SimpleFracture : public Modifier {
public:
  SimpleFracture(double alpha);

  virtual void initialize();
  virtual void evaluateStepOne(const int id_i, const int i);
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
  int m_indexBrokenNow;
  bool m_broken;
  arma::mat *m_data;
  std::unordered_map<int, int> *m_idToCol;
};
//------------------------------------------------------------------------------
}
#endif // SIMPLEFRACTURE_H
