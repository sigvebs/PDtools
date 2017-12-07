#ifndef CALCULATESTRESSSTRAIN_H
#define CALCULATESTRESSSTRAIN_H

#include "PDtools/CalculateProperties/calculateproperty.h"

namespace PDtools {
class Force;
//------------------------------------------------------------------------------
class CalculateStressStrain : public CalculateProperty {
public:
  CalculateStressStrain(vector<Force *> &forces, double E, double nu,
                        double delta, bool planeStress = false,
                        bool greenLagrange = false);

  virtual void initialize();

  virtual void update();

  void update2d();

  void update3d();

  void computeK(int id, int i);

  //       double weightFunction(const double dr0) const {return 1.;}
  double weightFunction(const double dr0) const { return m_delta / dr0; }
  //               double weightFunction(const double dr0) const {return
  //               m_delta*m_delta/(dr0*dr0);}
  //    double weightFunction(const double dr0) const {return
  //    exp(-2.*dr0/m_delta);}
  //       double weightFunction(const double dr0) const {return 1./(dr0*dr0);}
  //           double weightFunction(const double dr0) const {return (1. -
  //           1.01*dr0*dr0/(m_delta*m_delta));}
  //    double weightFunction(const double dr0) const {return 1.*(1. -
  //    1.01*dr0/m_delta);}
private:
  double m_E;
  double m_nu;
  double m_delta;
  double m_lambda;
  double m_mu;

  vector<Force *> m_forces;
  int m_indexStress[6];
  int m_indexStrain[6];
  int m_indexK[6];
  int m_nStressStrainElements;

  int m_indexConnected;
  int m_indexVolume;
  int m_indexVolumeScaling;
  int m_indexDr0;
  int m_indexBrokenNow;
  //    bool m_smallStrain = false;
  bool m_greenStrain = false;
  bool m_planeStress = false;
  bool m_analyticalK = false;
};
//------------------------------------------------------------------------------
}

#endif // CALCULATESTRESSSTRAIN_H
