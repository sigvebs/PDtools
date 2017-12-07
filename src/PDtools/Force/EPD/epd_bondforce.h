#ifndef EPD_BONDFORCE_H
#define EPD_BONDFORCE_H

#include "PDtools/Force/force.h"

namespace PDtools {
//------------------------------------------------------------------------------
class EPD_bondForce : public Force {
protected:
  int m_indexMicromodulus;
  int m_iOverlap;
  int m_iConnected;
  double m_dampCoeff;

  //    double weightFunction(const double dr0) const {return 1.;}
  //    double weightFunction(const double dr0) const {return 1.0/dr0;}
  double weightFunction(const double dr0) const { return m_delta / dr0; }
  //    double weightFunction(const double dr0) const {return
  //    exp(-3.*dr0/m_delta);}
  //    double weightFunction(const double dr0) const {return 1./(dr0*dr0);}
  //    double weightFunction(const double dr0) const {return m_delta - dr0;}
  //    double weightFunction(const double dr0) const {return (1. -
  //    1.001*dr0*dr0/(m_delta*m_delta));}
  //    double weightFunction(const double dr0) const {return 2.*(1. -
  //    1.01*dr0/m_delta);}
  void computeMicromodulus();

public:
  EPD_bondForce(PD_Particles &particles);

  virtual void calculateForces(const int id_i, const int i);

  virtual double calculateStableMass(const int id_a, const int a, double dt);

  virtual void initialize(double E, double nu, double delta, int dim, double h,
                          double lc);
};
//------------------------------------------------------------------------------
}
#endif // EPD_BONDFORCE_H
