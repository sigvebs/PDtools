#ifndef CALCULATESTRESSSTRAINEPD_H
#define CALCULATESTRESSSTRAINEPD_H

#include "PDtools/CalculateProperties/calculateproperty.h"
#include <vector>

namespace PDtools
{
class Force;
//------------------------------------------------------------------------------
class CalculateStressStrainEPD : public CalculateProperty
{
public:
    CalculateStressStrainEPD(vector<Force *> &forces,
                             double E, double nu, double delta,
                             bool planeStress=false,
                             bool greenLagrange=false);


    virtual void
    initialize();

    virtual void
    clean();

    virtual void
    update();

    void
    computeK(int id, int i);

    //       double weightFunction(const double dr0) const {return 1.;}
    double weightFunction(const double dr0) const {return m_delta/dr0;}
//    double weightFunction(const double dr0) const {return exp(-3.*dr0/m_delta);}
    //              double weightFunction(const double dr0) const {return 5.*(1. - 1.01*dr0/m_delta);}
    //       double weightFunction(const double dr0) const {return 1./(dr0*dr0);}
    //              double weightFunction(const double dr0) const {return 2.*(1. - 1.01*dr0/m_delta);}
    //       double weightFunction(const double dr0) const {return (1. - 1.001*dr0*dr0/(m_delta*m_delta));}
    //              double weightFunction(const double dr0) const {return 1.*(1. - 1.001*dr0*dr0/(m_delta*m_delta));}
    //       double weightFunction(const double dr0) const {return exp(-5.*dr0/m_delta);}
//    double weightFunction(const double dr0) const {return 0.5*(1. - 1.01*dr0/m_delta);}
private:
    vector<Force *> m_forces;
    int m_indexStress[6];
    int m_indexStrain[6];
    int m_indexK[6];
    int m_nStressStrainElements;

    int m_indexConnected;
    int m_indexVolume;
    int m_indexVolumeScaling;
    int m_indexDr0;
    int m_indexBrokenNow;
    int m_iOverlap;
    bool m_smallStrain = false;
    bool m_greenStrain = false;
    bool m_planeStress = false;

    double m_E;
    double m_nu;
    double m_delta;
};
//------------------------------------------------------------------------------
}


#endif // CALCULATESTRESSSTRAINEPD_H
