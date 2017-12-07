#include "calculateproperty.h"

namespace PDtools {
//------------------------------------------------------------------------------
CalculateProperty::CalculateProperty(string _type) : type(_type) {}
//------------------------------------------------------------------------------
CalculateProperty::~CalculateProperty() {}
//------------------------------------------------------------------------------
void CalculateProperty::initialize() {}
//------------------------------------------------------------------------------
int CalculateProperty::updateFrequency() const { return m_updateFrequency; }
//------------------------------------------------------------------------------
void CalculateProperty::setUpdateFrquency(int updateFrquency) {
  m_updateFrequency = updateFrquency;
}
//------------------------------------------------------------------------------
void CalculateProperty::clean() {}
//------------------------------------------------------------------------------
int CalculateProperty::dim() const { return m_dim; }
//------------------------------------------------------------------------------
void CalculateProperty::setDim(int dim) { m_dim = dim; }
//------------------------------------------------------------------------------
void CalculateProperty::setParticles(PD_Particles &particles) {
  m_particles = &particles;
}
//------------------------------------------------------------------------------
}
