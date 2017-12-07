#ifndef STRAINBOUNDARY_H
#define STRAINBOUNDARY_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class StrainBoundary : public Modifier {
public:
  StrainBoundary(double strainrate, double strainOrientation, double nu,
                 pair<double, double> boundary, int boundaryOrientation,
                 double dt, bool isStatic);

  virtual void registerParticleParameters();
  virtual void evaluateStepOne();
  virtual void initialize();
  virtual void staticEvaluation();

private:
  double m_strainrate;
  int m_strainOritentation;
  pair<double, double> m_boundary;
  int m_boundaryOrientation;
  vector<int> m_otherDirections;
  double m_dt;
  double m_time;
  double m_nu;
  bool m_isStatic;
  bool m_usingUnbreakableBorder;
};
//------------------------------------------------------------------------------
}
#endif // STRAINBOUNDARY_H
