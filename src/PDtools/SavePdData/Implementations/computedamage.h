#ifndef COMPUTEDAMAGE_H
#define COMPUTEDAMAGE_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools {

class ComputeDamage : public ComputeProperty {
public:
  ComputeDamage(PD_Particles &particles);
  ~ComputeDamage();

  virtual void update(const int id_i, const int i);
  virtual void init(const int id_i, const int i);

private:
  int m_indexDamage;
  int m_indexMaxPdConnections;
  int m_indexConnected;
  arma::mat *m_data;
};
//------------------------------------------------------------------------------
}
#endif // COMPUTEDAMAGE_H
