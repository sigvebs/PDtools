#ifndef SOLVER_H
#define SOLVER_H

struct PD_Particle;
class Domain;

using namespace std;

namespace PDtools
{
//------------------------------------------------------------------------------
class Solver
{
public:
    Solver(Domain & domain);
    virtual void solve();
    virtual void stepForward( int i, double t ) = 0;
    void setSingleParticleForces(){;}
    void setSingleParticleModifiers(){;}
protected:
    Domain & domain;
    double dt;
    double t_start;
    int steps;

    //    vector<modifier*> singleParticeModifieres;
    //    vector<force*> singleParticeForces;
};
//------------------------------------------------------------------------------
}
#endif // SOLVER_H
