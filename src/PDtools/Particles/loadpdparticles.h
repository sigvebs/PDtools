#ifndef LOADPDPARTICLES_H
#define LOADPDPARTICLES_H

#include "loadparticles.h"

namespace PDtools {
class PD_Particles;

//------------------------------------------------------------------------------
class LoadPdParticles : public LoadParticles {
public:
  LoadPdParticles();

  PD_Particles load(string loadPath, string format, bool bin = false,
                    unordered_map<string, int> loadParameters = {{}});

protected:
  const vector<string> m_PdBasicParameters = {"id",  "x",   "y",  "z",
                                              "v_x", "v_y", "v_z"};
  virtual void loadBody(PD_Particles &particles, std::fstream &rawData,
                        unordered_map<string, int> parameters);
  virtual void loadBinaryBody(PD_Particles &particles, FILE *rawData,
                              unordered_map<string, int> parameters);
};
//------------------------------------------------------------------------------
PD_Particles load_pd(string loadPath);
PD_Particles load_pd(string loadPath, Grid &grid);
PD_Particles load_pd(string loadPath,
                     unordered_map<string, double> particleParameters);
//------------------------------------------------------------------------------
}
#endif // LOADPDPARTICLES_H
