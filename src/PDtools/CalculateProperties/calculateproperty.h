#ifndef CALCULATEPROERTY_H
#define CALCULATEPROERTY_H

#include<armadillo>
using namespace std;

namespace PDtools
{
class PD_Particles;

//------------------------------------------------------------------------------
class CalculateProperty
{
public:
    const string type;

    CalculateProperty(string _type);
    ~CalculateProperty();

    virtual void
    initialize();

    int
    updateFrquency() const;

    void
    setUpdateFrquency(int updateFrquency);

    virtual void
    clean();

    virtual void
    update() = 0;

    int
    dim() const;

    void
    setDim(int dim);

    void setParticles(PD_Particles &particles);

protected:
    PD_Particles *m_particles;
    int m_dim;
    int m_updateFrquency = 1;
};
//------------------------------------------------------------------------------
}
#endif // CALCULATEPROERTY_H