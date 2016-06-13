#ifndef MOHRCOULOMBMAXFRACTUREWEIGHTEDADR_H
#define MOHRCOULOMBMAXFRACTUREWEIGHTEDADR_H

#include "PDtools/Modfiers/Implementation/FractureCriterion/mohrcoulombmaxfractureweighted.h"

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
class MohrCoulombMaxFractureWeightedAdr : public MohrCoulombMaxFractureWeighted
{
public:
    MohrCoulombMaxFractureWeightedAdr(double mu, double C, double T, double Wc, double Bc);

    virtual void
    evaluateStepOne();

    virtual void
    evaluateStepOne(const int id_i, const int i);


    virtual void
    evaluateStepTwo(const int id_i, const int i);

    virtual void
    evaluateStepTwoPost();
};
//------------------------------------------------------------------------------
}
#endif // MOHRCOULOMBMAXFRACTUREWEIGHTEDADR_H
