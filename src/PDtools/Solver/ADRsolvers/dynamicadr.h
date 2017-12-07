#ifndef DYNAMICADR_H
#define DYNAMICADR_H

#include "PDtools/Solver/adr.h"

namespace PDtools {
//------------------------------------------------------------------------------
class dynamicADR : public ADR {
public:
  dynamicADR();
  virtual void solve();
  virtual void stepForward(int timeStep);
};
//------------------------------------------------------------------------------
}

#endif // DYNAMICADR_H
