#ifndef CALCULATESTRESS_H
#define CALCULATESTRESS_H

#include "PDtools/CalculateProperties/calculateproperty.h"

namespace PDtools {
class Force;
//------------------------------------------------------------------------------

class CalculateStress : public CalculateProperty {
public:
  CalculateStress(vector<Force *> &forces);

  virtual void initialize();

  virtual void clean();

  virtual void update();

private:
  vector<Force *> m_forces;
  int m_indexStress[6];
  int m_nStressElements;
};
//------------------------------------------------------------------------------
}
#endif // CALCULATESTRESS_H
