#ifndef MOVEPARTICLES_H
#define MOVEPARTICLES_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools {
//------------------------------------------------------------------------------
class MoveParticles : public Modifier {
public:
  MoveParticles(double velAmplitude, double velOrientation,
                pair<double, double> boundary, int boundaryOrientation,
                double dt, bool isStatic);

  virtual void registerParticleParameters();
  virtual void evaluateStepOne();
  virtual void initialize();
  virtual void staticEvaluation();

private:
  double m_velAmplitude;
  int m_velOritentation;
  pair<double, double> m_boundary;
  int m_boundaryOrientation;
  double m_dt;
  double m_time;
  bool m_isStatic;
  bool m_usingUnbreakableBorder;
};
//------------------------------------------------------------------------------
}

#endif // MOVEPARTICLES_H
