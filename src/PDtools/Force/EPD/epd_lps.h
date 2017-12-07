#ifndef EPD_LPS_H
#define EPD_LPS_H

#include "PDtools/Force/force.h"

namespace PDtools {
//------------------------------------------------------------------------------
class EPD_LPS : public Force {
protected:
  bool m_planeStress;
  int m_dim = 3;
  int m_iMicromodulus;
  int m_iTheta;
  int m_iThetaNew;
  int m_iVolume;
  int m_iDr0;
  int m_iOverlap;
  int m_iStretch;
  int m_iForceScalingDilation;
  int m_iForceScalingBond;
  int m_iConnected;
  int m_iMass;
  int m_iCompute;
  int m_indexBrokenNow;

  double m_k;
  double m_mu;
  double m_delta;
  double m_nu;
  double m_alpha;
  double m_c;
  double m_t;

public:
  EPD_LPS(PD_Particles &particles, bool planeStress = false);

  virtual void calculateForces(const int id, const int i);

  virtual double calculatePotentialEnergyDensity(const int id_i, const int i);

  double computeDilation(const int id_i, const int i);

  virtual void calculatePotentialEnergy(const int id_i, const int i,
                                        int indexPotential);
  virtual void calculateStress(const int id_i, const int i,
                               const int (&indexStress)[6]);

  virtual void updateState(int id, int i);

  virtual double calculateStableMass(const int id_a, const int a, double dt);
  virtual void initialize(double E, double nu, double delta, int dim, double h,
                          double lc);

  void calculateWeightedVolume();
};
//------------------------------------------------------------------------------
}
#endif // EPD_LPS_H
