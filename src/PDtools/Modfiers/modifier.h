#ifndef MODIFIER_H
#define MODIFIER_H

#include "config.h"

namespace PDtools {
class PD_Particles;
class Grid;

//------------------------------------------------------------------------------
class Modifier {
protected:
  PD_Particles *m_particles;
  Grid *m_grid;
  bool m_state = false;
  vector<int> m_localParticleIds;
  std::vector<std::string> m_initialGhostParameters;
  std::vector<std::string> m_ghostParameters;
  std::vector<pair<string, int>>
      m_neededProperties; // name and update frequency
  int m_dim;
  int m_myRank = 0;
  int m_nCores = 1;

  bool m_hasStepOne = false;
  bool m_hasUpdateOne = false;
  bool m_hasStepTwo = false;

public:
  Modifier();
  virtual ~Modifier();

  virtual void registerParticleParameters();
  virtual void evaluateStepOne(const int id, const int i);
  virtual void updateStepOne(const int id, const int i);
  virtual void evaluateStepTwo(const int id, const int i);
  virtual void evaluateStepOne();
  virtual void evaluateStepOnePost();
  virtual void evaluateStepTwo();
  virtual void staticEvaluation();
  virtual void initialize();
  void setParticles(PD_Particles &particles);
  bool state();

  const std::vector<std::string> &initalGhostDependencies();

  const std::vector<std::string> &ghostDependencies();

  void addToList(int id);

  bool removeFromList(const int id);

  std::vector<std::pair<std::string, int>> neededProperties() const;
  int dim() const;
  void setDim(int dim);

  void setGrid(Grid *grid);
  bool hasStepOne() const;
  bool hasStepTwo() const;
  bool hasUpdateOne() const;
};
//------------------------------------------------------------------------------
}
#endif // MODIFIER_H
