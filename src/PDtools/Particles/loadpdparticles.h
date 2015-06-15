#ifndef LOADPDPARTICLES_H
#define LOADPDPARTICLES_H

#include "loadparticles.h"

using namespace std;

namespace PDtools
{
// Forward declerations
class PD_Particles;

//------------------------------------------------------------------------------
class LoadPdParticles: public T_LoadParticles<PD_Particles>
{
public:
    LoadPdParticles();
protected:
    const vector<string> m_PdBasicParameters = {"id", "x", "y", "z",
                                                "v_x", "v_y", "v_z"};

    virtual void loadBody(PD_Particles &particles, fstream &rawData,
                  unordered_map <string, int> parameters);

    virtual void loadBinaryBody(PD_Particles &particles, FILE *rawData,
                        unordered_map <string, int> parameters);
};
//------------------------------------------------------------------------------
PD_Particles load_pd(string loadPath);
PD_Particles load_pd(string loadPath,
                     unordered_map<string, double> particleParameters);
//------------------------------------------------------------------------------
}
#endif // LOADPDPARTICLES_H
