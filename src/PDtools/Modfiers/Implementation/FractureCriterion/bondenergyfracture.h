#ifndef BONDENERGYFRACTURE_H
#define BONDENERGYFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
class Force;
//------------------------------------------------------------------------------
class BondEnergyFracture : public Modifier {
public:
  BondEnergyFracture(double delta, double G0, vector<Force *> &forces,
                     double h = 1.);

  virtual void initialize();
  virtual void evaluateStepOne(const int id_i, const int i);
  virtual void evaluateStepTwo();

protected:
  vector<Force *> &m_forces;
  double m_delta;
  double m_wc; // Critical energy density
  int m_indexStretch;
  int m_indexUnbreakable;
  int m_indexConnected;
  int m_indexDr0;
  int m_broken;
  int m_indexForceScaling;
  int m_indexMicromodulus;
  int m_indexBrokenNow;
  double m_G;
  double m_h;
  arma::mat *m_data;
  ivec *m_idToCol;
};
//------------------------------------------------------------------------------
}
#endif // BONDENERGYFRACTURE_H
