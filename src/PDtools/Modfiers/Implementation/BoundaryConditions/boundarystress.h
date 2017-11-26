#ifndef BOUNDARYSTRESS_H
#define BOUNDARYSTRESS_H

#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class BoundaryStress : public Modifier
{
public:
    BoundaryStress(double appliedStress, double stressOrientation,
                   pair<double, double> boundary, int boundaryOrientation,
                   int steps, double delta);

  virtual void evaluateStepOne();
  virtual void evaluateStepTwo();
  virtual void initialize();
  virtual void staticEvaluation();

private:
  double m_appliedStress;
  int m_stressOritentation;
  pair<double, double> m_boundary;
  int m_boundaryOrientation;
  vector<int> m_otherAxis;
  int m_indexVolume ;
  int m_steps;
  int m_indexRadius;
  double m_incrementalStress;
  int m_inceremental;
  double m_delta;
};
//------------------------------------------------------------------------------
}
#endif // BOUNDARYSTRESS_H
