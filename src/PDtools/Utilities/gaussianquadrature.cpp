#include "gaussianquadrature.h"

namespace PDtools
{
//------------------------------------------------------------------------------
GaussLegendreQuad::GaussLegendreQuad(int _dim, int _n):
    dim(_dim), n(_n)
{
    m_shapeFunction = {0, 0, 0, 0};
    m_x.reserve(_n);
    m_w.reserve(_n);

    // For quads
    switch (_n) {
    case 1:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x1[i]);
            m_w.push_back(w1[i]);
        }
        break;
    case 2:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x2[i]);
            m_w.push_back(w2[i]);
        }
        break;
    case 3:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x3[i]);
            m_w.push_back(w3[i]);
        }
        break;
    case 4:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x4[i]);
            m_w.push_back(w4[i]);
        }
        break;
    case 5:
        for(int i=0; i<_n; i++) {
            m_x.push_back(x5[i]);
            m_w.push_back(w5[i]);
        }
        break;
    default:
        cerr << "Degree of Gauss-Legendre polynomial not implemented: " << _n << endl;
        exit(EXIT_FAILURE);
        break;
    }

    //--------------------------------------------------------------------------
    // For 2d
    //--------------------------------------------------------------------------
    const int shapeFunction_size = 4;
    m_nPoints = n*n;
    m_gaussianPoints_2d = zeros(n*n, dim);
    m_gaussianShapFunction_2d = zeros(n*n, shapeFunction_size);
    m_weights_2d = zeros(n*n);

    int counter = 0;
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            const double x = m_x[i];
            const double y = m_x[j];
            m_gaussianPoints_2d(counter, 0) = x;
            m_gaussianPoints_2d(counter, 1) = y;
            m_weights_2d(counter) = m_w[i]*m_w[j];
            vector<double> shape = shapeFunction(x, y);
            for(int k=0; k<shape.size(); k++) {
                m_gaussianShapFunction_2d(counter, k) = shape[k];
            }

            // Helper arrays for calculating the Jakobi determinant
            vector<double> a1 = {-(1-m_x[i]),   1-m_x[i] , 1 + m_x[i], -(1 + m_x[i])};
            vector<double> a2 = {-(1-m_x[j]), -(1+m_x[j]), 1 + m_x[j],  (1 - m_x[j])};
            for(int k=0; k<4; k++) {
                a1[k] *= 1./4.;
                a2[k] *= 1./4.;
            }
            m_jakobi_1.push_back(a1);
            m_jakobi_2.push_back(a2);

            counter++;
        }
    }

//    for(int i=0; i<_n; i++) {
//        vector<double> a1 = {-(1-m_w[i]),   1-m_w[i] , 1 + m_w[i], -(1 + m_w[i])};
//        vector<double> a2 = {-(1-m_w[i]), -(1+m_w[i]), 1 + m_w[i],  (1 - m_w[i])};
//        for(int j=0; j<4; j++) {
//            a1[j] *= 1./4.;
//            a2[j] *= 1./4.;
//        }
//        m_jakobi_1.push_back(a1);
//        m_jakobi_2.push_back(a2);
//    }

}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::shapeFunction(double x, double y)
{
    switch(dim) {
    case 2:
        m_shapeFunction[0] = 1./4.*(1-x)*(1-y);
        m_shapeFunction[1] = 1./4.*(1+x)*(1-y);
        m_shapeFunction[2] = 1./4.*(1+x)*(1+y);
        m_shapeFunction[3] = 1./4.*(1-x)*(1+y);
        break;
    }

    return m_shapeFunction;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::x() const
{
    return m_x;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreQuad::w() const
{
    return m_w;
}
//------------------------------------------------------------------------------
mat GaussLegendreQuad::gaussianPoints_2d()
{
    return m_gaussianPoints_2d;
}
//------------------------------------------------------------------------------
vec GaussLegendreQuad::gaussianWeights_2d()
{
    return m_weights_2d;
}
//------------------------------------------------------------------------------
mat GaussLegendreQuad::shapeFunction_2d()
{
    return m_gaussianShapFunction_2d;
}

int GaussLegendreQuad::nPoints() const
{
    return m_nPoints;
}
//------------------------------------------------------------------------------
vector<vector<double> > GaussLegendreQuad::jakobi_1() const
{
    return m_jakobi_1;
}
//------------------------------------------------------------------------------
vector<vector<double> > GaussLegendreQuad::jakobi_2() const
{
    return m_jakobi_2;
}
//------------------------------------------------------------------------------
GaussLegendreTriangle::GaussLegendreTriangle(int _dim, int _n):
    dim(_dim), n(_n)
{
    m_shapeFunction = {0, 0, 0};
    m_x.reserve(_n);
    m_w.reserve(_n);


    // For triangles
    switch (_n) {
    case 1:
        m_nTriPoints = 1;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x1[i]);
            m_y.push_back(y1[i]);
            m_w.push_back(w1[i]);
        }
        break;
    case 2:
        m_nTriPoints = 3;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x2[i]);
            m_y.push_back(y2[i]);
            m_w.push_back(w2[i]);
        }
        break;
    case 3:
        m_nTriPoints = 4;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x3[i]);
            m_y.push_back(y3[i]);
            m_w.push_back(w3[i]);
        }
        break;
    case 4:
        m_nTriPoints = 6;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x4[i]);
            m_y.push_back(y4[i]);
            m_w.push_back(w4[i]);
        }
        break;
    case 5:
        m_nTriPoints = 7;
        for(int i=0; i<m_nTriPoints; i++) {
            m_x.push_back(x5[i]);
            m_y.push_back(y5[i]);
            m_w.push_back(w5[i]);
        }
        break;
    default:
        cerr << "Degree of Gauss-Legendre polynomial for triangles not implemented for degree : " << _n << endl;
        exit(EXIT_FAILURE);
        break;
    }
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::shapeFunction(double x, double y)
{
    switch(dim) {
    case 2:
        m_shapeFunction[0] = 1. - x - y;
        m_shapeFunction[1] = x;
        m_shapeFunction[2] = y;
        break;
    }

    return m_shapeFunction;
}
//------------------------------------------------------------------------------
int GaussLegendreTriangle::nTriPoints() const
{
    return m_nTriPoints;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::x() const
{
    return m_x;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::y() const
{
    return m_y;
}
//------------------------------------------------------------------------------
vector<double> GaussLegendreTriangle::w() const
{
    return m_w;
}
//------------------------------------------------------------------------------
orderedGrid::orderedGrid(int _dim, int _n):
    dim(_dim), n(_n)
{
    m_shapeFunction = {0, 0, 0, 0};
    m_x.reserve(_n);
    m_w.reserve(_n);

    const double dx = 2./_n;
    for(int i=0; i<_n; i++) {
        double x = i*dx + 0.5*dx - 1.;
        m_x.push_back(x);
        m_w.push_back(dx);
    }

    //--------------------------------------------------------------------------
    // For 2d
    //--------------------------------------------------------------------------
    const int shapeFunction_size = 4;
    m_nPoints = n*n;
    m_gaussianPoints_2d = zeros(n*n, dim);
    m_gaussianShapFunction_2d = zeros(n*n, shapeFunction_size);
    m_weights_2d = zeros(n*n);

    int counter = 0;
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            const double x = m_x[i];
            const double y = m_x[j];
            m_gaussianPoints_2d(counter, 0) = x;
            m_gaussianPoints_2d(counter, 1) = y;
            m_weights_2d(counter) = m_w[i]*m_w[j];
            vector<double> shape = shapeFunction(x, y);
            for(int k=0; k<shape.size(); k++) {
                m_gaussianShapFunction_2d(counter, k) = shape[k];
            }

            // Helper arrays for calculating the Jakobi determinant
            vector<double> a1 = {-(1-m_x[i]),   1-m_x[i] , 1 + m_x[i], -(1 + m_x[i])};
            vector<double> a2 = {-(1-m_x[j]), -(1+m_x[j]), 1 + m_x[j],  (1 - m_x[j])};
            for(int k=0; k<4; k++) {
                a1[k] *= 1./4.;
                a2[k] *= 1./4.;
            }
            m_jakobi_1.push_back(a1);
            m_jakobi_2.push_back(a2);

            counter++;
        }
    }
}
//------------------------------------------------------------------------------
vector<double> orderedGrid::shapeFunction(double x, double y)
{
    switch(dim) {
    case 2:
        m_shapeFunction[0] = 1./4.*(1-x)*(1-y);
        m_shapeFunction[1] = 1./4.*(1+x)*(1-y);
        m_shapeFunction[2] = 1./4.*(1+x)*(1+y);
        m_shapeFunction[3] = 1./4.*(1-x)*(1+y);
        break;
    }

    return m_shapeFunction;
}
//------------------------------------------------------------------------------
vector<double> orderedGrid::x() const
{
    return m_x;
}
//------------------------------------------------------------------------------
vector<double> orderedGrid::w() const
{
    return m_w;
}
//------------------------------------------------------------------------------
mat orderedGrid::gaussianPoints_2d()
{
    return m_gaussianPoints_2d;
}
//------------------------------------------------------------------------------
vec orderedGrid::gaussianWeights_2d()
{
    return m_weights_2d;
}
//------------------------------------------------------------------------------
mat orderedGrid::shapeFunction_2d()
{
    return m_gaussianShapFunction_2d;
}
//------------------------------------------------------------------------------
int orderedGrid::nPoints() const
{
    return m_nPoints;
}
//------------------------------------------------------------------------------
vector<vector<double> > orderedGrid::jakobi_1() const
{
    return m_jakobi_1;
}
//------------------------------------------------------------------------------
vector<vector<double> > orderedGrid::jakobi_2() const
{
    return m_jakobi_2;
}
//------------------------------------------------------------------------------
}
