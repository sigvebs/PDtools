#ifndef MOHRCOULOMBMAXFRACTURE_H
#define MOHRCOULOMBMAXFRACTURE_H

#include "PDtools/Modfiers/modifier.h"

//------------------------------------------------------------------------------
namespace PDtools {
//------------------------------------------------------------------------------
class MohrCoulombMaxFracture : public Modifier {
public:
  MohrCoulombMaxFracture(double mu, double C, double T);

  virtual void registerParticleParameters();
  virtual void initialize();
  virtual void evaluateStepOne(const int id_i, const int i);
  virtual void evaluateStepOnePost();
  virtual void evaluateStepTwo(const int id_i, const int i);
  virtual void evaluateStepTwo();

private:
  double m_C;
  double m_T;
  double m_d;
  double m_phi;
  arma::mat *m_data;
  std::unordered_map<int, int> *m_idToCol;
  unordered_map<int, vector<int>> m_brokenParticles;

  int m_indexStress[6];
  int m_indexNormal[3];
  int m_indexUnbreakable;
  int m_indexConnected;
  int m_indexCompute;
  int m_indexStressCenter;
  int m_indexBrokenNow;
  int m_indexBroken;
  int m_indexRadius;
  int m_indexDamage;
  bool m_broken;
  double m_cosTheta;
  double m_sinTheta;

  void exchangeBrokenParticlesMPI();
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBMAXFRACTURE_H
