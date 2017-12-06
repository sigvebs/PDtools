#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <armadillo>

using namespace std;
namespace PDtools
{

//------------------------------------------------------------------------------
class Domain
{
public:
    const int dim;

    Domain(int dim, vector <pair<double, double>> boundaries);
    ~Domain() {}

    void setBoundaries(vector <pair<double, double>> boundaries);

    const vector <pair<double, double>> & boundaries() {
        return m_boundaries;
    }

    void periodicBoundaries(arma::ivec3 pb) {
        m_periodicBoundaries = pb;
    }

    const arma::ivec3 & periodicBoundaries() const {
        return m_periodicBoundaries;
    }

    const vector <double> boundaryLength() const {
        return m_boundaryLength;
    }

    const vector <pair<double, double>> boundaries() const {
        return m_boundaries;
    }

private:
    arma::ivec3 m_periodicBoundaries;
    vector <pair<double, double>> m_boundaries;
    vector <double> m_boundaryLength;
};
//------------------------------------------------------------------------------

}
#endif // DOMAIN_H
