#ifndef COMPUTEMAXSTRETCH_H
#define COMPUTEMAXSTRETCH_H

#include "PDtools/SavePdData/savepddata.h"
namespace PDtools {

//------------------------------------------------------------------------------
class ComputeMaxStretch : public ComputeProperty {
public:
  ComputeMaxStretch(PD_Particles &particles);
  ~ComputeMaxStretch();

  virtual void update(const int id_i, const int i);
  virtual void init(const int id_i, const int i);

private:
  int m_indexDamage;
  int m_indexStretch;
  int m_indexMaxStretch;
  arma::mat *m_data;
};
//------------------------------------------------------------------------------
}
#endif // COMPUTEMAXSTRETCH_H
