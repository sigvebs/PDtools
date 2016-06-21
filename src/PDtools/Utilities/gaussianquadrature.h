#ifndef GAUSSIANQUADRATURE_H
#define GAUSSIANQUADRATURE_H

#include <iostream>
#include <vector>
#include <armadillo>

namespace PDtools
{
using namespace arma;
using namespace std;
//------------------------------------------------------------------------------
class GaussLegendreQuad
{
public:
    GaussLegendreQuad(int _dim, int _n);
    const int dim;
    const int n;

    vector<double> shapeFunction(double x, double y);

    vector<double> x() const;
    vector<double> w() const;

    mat gaussianPoints_2d();
    vec gaussianWeights_2d();
    mat shapeFunction_2d();

    int nPoints() const;

protected:
    vector<double> m_shapeFunction;
    vector<double> m_x;
    vector<double> m_w;
    int m_nPoints;
    mat m_gaussianPoints_2d;
    mat m_gaussianShapFunction_2d;
    vec m_weights_2d;

    // For quads
    const double x1[1] = {0.};
    const double w1[1] = {2.0000000000000000000000000};

    const double x2[2] = {0.5773502691896257645091488, -0.5773502691896257645091488};
    const double w2[2] = {1.0000000000000000000000000, 1.0000000000000000000000000};

    const double x3[3] = {0.0000000000000000000000000, 0.7745966692414833770358531, -0.7745966692414833770358531};
    const double w3[3] = {0.8888888888888888888888889, 0.5555555555555555555555556, 0.5555555555555555555555556};

    const double x4[4] = {0.3399810435848562648026658, 0.8611363115940525752239465, -0.3399810435848562648026658, -0.8611363115940525752239465};
    const double w4[4] = {0.6521451548625461426269361, 0.3478548451374538573730639, 0.6521451548625461426269361, 0.3478548451374538573730639};

    const double x5[5] = {0.0000000000000000000000000, 0.5384693101056830910363144, 0.9061798459386639927976269, -0.5384693101056830910363144, -0.9061798459386639927976269};
    const double w5[5] = {0.5688888888888888888888889, 0.4786286704993664680412915, 0.2369268850561890875142640, 0.4786286704993664680412915, 0.2369268850561890875142640};
};
//------------------------------------------------------------------------------
class GaussLegendreTriangle
{
public:
    GaussLegendreTriangle(int _dim, int _n);
    const int dim;
    const int n;

    vector<double> shapeFunction(double x, double y);

    int nTriPoints() const;

    vector<double> x() const;
    vector<double> y() const;
    vector<double> w() const;

protected:
    vector<double> m_shapeFunction;
    vector<double> m_x;
    vector<double> m_y;
    vector<double> m_w;
    int m_nTriPoints;

    // For triangles
    const double x1[1] = {1./3.};
    const double y1[1] = {1./3.};
    const double w1[1] = {1.};

    const double x2[3] = {1./6., 2./3., 1./6.};
    const double y2[3] = {1./6., 1./6., 2./3.};
    const double w2[3] = {1./3, 1./3, 1./3};

    const double x3[4] = {1./3., 1./5., 1./5., 3./5.};
    const double y3[4] = {1./3., 3./5., 1./5., 1./5.};
    const double w3[4] = {-27./48., 25./48., 25./48., 25./48.};

    const double x4[6] = {0.44594849091597, 0.44594849091597, 0.10810301816807, 0.09157621350977, 0.09157621350977, 0.81684757298046};
    const double y4[6] = {0.44594849091597, 0.10810301816807, 0.44594849091597, 0.09157621350977, 0.81684757298046, 0.09157621350977};
    const double w4[6] = {0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532};

    const double x5[7] = {0.33333333333333, 0.47014206410511, 0.47014206410511, 0.05971587178977, 0.10128650732346, 0.10128650732346, 0.79742698535309};
    const double y5[7] = {0.33333333333333, 0.47014206410511, 0.05971587178977, 0.47014206410511, 0.10128650732346, 0.79742698535309, 0.10128650732346};
    const double w5[7] = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
};
//------------------------------------------------------------------------------
}
#endif // GAUSSIANQUADRATURE_H
