#ifndef VONMISESFRACTURE_H
#define VONMISESFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

//------------------------------------------------------------------------------
namespace PDtools {
//------------------------------------------------------------------------------
class VonMisesFracture : public Modifier {
public:
  VonMisesFracture(double sigma_y);

  virtual void registerParticleParameters();

  virtual void evaluateStepOne(const int id_i, const int i);

protected:
  double m_sigma_y;
  double m_sigma_y2;
  arma::mat *m_data;
  ivec *m_idToCol;

  int m_indexUnbreakable;
  int m_indexConnected;
  int m_indexCompute;
  int m_indexBrokenNow;
  int m_indexStress[6];
};
//------------------------------------------------------------------------------
}

#endif // VONMISESFRACTURE_H
