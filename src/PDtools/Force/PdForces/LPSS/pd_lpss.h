#ifndef PD_LPSS_H
#define PD_LPSS_H

#include "PDtools/Force/force.h"

namespace PDtools {

//------------------------------------------------------------------------------
class PD_LPSS : public Force {
protected:
  bool m_planeStress;
  int m_iMicromodulus;
  int m_iTheta;
  int m_iThetaNew;
  int m_iVolume;
  int m_iDr0;
  int m_iVolumeScaling;
  int m_iStretch;
  int m_iForceScalingDilation;
  int m_iForceScalingBond;
  int m_iConnected;
  int m_iMass;
  int m_iCompute;
  int m_iBrokenNow;
  bool m_analyticalM = false;
  double m_k;
  double m_mu;

  // Pointers
  double *m_mass;
  double *m_theta;
  double *m_volume;
  double *m_x;
  double *m_y;
  double *m_z;
  double *m_x0;
  double *m_y0;
  double *m_z0;
  double *m_Fx;
  double *m_Fy;
  double *m_Fz;

  int m_iK[6]; // Shape matrix
  int m_iR[9]; // Rotation matrix
  int m_iStress[6];

  double weightFunction(const double dr0) const { return 1. / (dr0 * dr0); }
  //    double weightFunction(const double dr0) const {return m_delta/dr0;}
public:
  PD_LPSS(PD_Particles &particles, bool planeStress = false);

  virtual void calculateForces(const int id, const int i);

  virtual double calculatePotentialEnergyDensity(const int id_i, const int i);

  virtual void calculatePotentialEnergy(const int id_i, const int i,
                                        int indexPotential);

  virtual void evaluateStepOne();

  virtual void evaluateStepOne(const int id_i, const int i);

  virtual double calculateStableMass(const int id_a, const int a, double dt);

  void calculateWeightedVolume();

  void computeMandK(int id_i, int i);

  virtual void initialize(double E, double nu, double delta, int dim, double h,
                          double lc);

  double computeDilation(const int id_i, const int i);

  void updateWeightedVolume(int id_i, int i);
};
//------------------------------------------------------------------------------
}
#endif // PD_LPSS_H
