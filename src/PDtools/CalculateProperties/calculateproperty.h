#ifndef CALCULATEPROERTY_H
#define CALCULATEPROERTY_H

#include "config.h"

namespace PDtools {
class PD_Particles;

//------------------------------------------------------------------------------
class CalculateProperty {
public:
  const string type;

  CalculateProperty(string _type);
  virtual ~CalculateProperty();

  virtual void initialize();
  int updateFrequency() const;
  void setUpdateFrquency(int updateFrequency);
  virtual void clean();
  virtual void update() = 0;
  int dim() const;
  void setDim(int dim);
  void setParticles(PD_Particles &particles);

protected:
  PD_Particles *m_particles;
  int m_dim;
  int m_updateFrequency = 1;
};
//------------------------------------------------------------------------------
}
#endif // CALCULATEPROERTY_H
