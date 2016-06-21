#ifndef PD_ELEMENT_H
#define PD_ELEMENT_H

#include <vector>
#include <array>
#include <armadillo>

//#include "boostgeometry_settings.h"

//------------------------------------------------------------------------------
namespace PDtools
{
using namespace std;

class PdMesh;
class PD_Particles;
class Grid;

//enum elementTYpes {TRIANGLE=3, QUAD=4};

//------------------------------------------------------------------------------
template <size_t T>
class PD_element
{
public:
    PD_element() {}
    PD_element(size_t _id, array<size_t, T> _verticeIds): m_id(_id), m_verticeIds(_verticeIds) {}
    const array<size_t, T> & verticeIds() const {return m_verticeIds;}

    size_t id() const {return m_id;}

    arma::mat guassianQuadraturePoints_initial() const {return m_guassianQuadraturePoints_initial;}
    void setGuassianQuadraturePoints_initial(const arma::mat &guassianQuadraturePoints_initial)
    {m_guassianQuadraturePoints_initial = guassianQuadraturePoints_initial;}

    arma::mat guassianQuadraturePoints() const {return m_guassianQuadraturePoints;}
    void setGuassianQuadraturePoints(const arma::mat &guassianQuadraturePoints)
    {m_guassianQuadraturePoints = guassianQuadraturePoints;}

    arma::mat guassianQuadratureWeights() const {return m_guassianQuadratureWeights;}
    void setGuassianQuadratureWeights(const arma::mat &guassianQuadratureWeights)
    {m_guassianQuadratureWeights = guassianQuadratureWeights;}

protected:
    const size_t m_id;
    array<size_t, T> m_verticeIds;
    arma::mat m_guassianQuadraturePoints_initial;
    arma::mat m_guassianQuadraturePoints;
    arma::vec m_guassianQuadratureWeights;
};
typedef PD_element<2> PD_lineElement;
typedef PD_element<3> PD_triElement;
typedef PD_element<4> PD_quadElement;

PD_Particles initializeElementPd(const PdMesh &msh, const Grid & grid, const size_t quadratureDegree);
//------------------------------------------------------------------------------
//class PD_elementsAndNodes
//{
//public:
//    PD_elementsAndNodes(Grid & grid);

//    void initialize(const PdMesh &msh);

//protected:
//    const int m_dim;
//    const Grid &m_grid;
//    PD_Particles m_particles;
//    vector<PD_triElement> m_triElements;
//    vector<PD_quadElement> m_quadElements;
//};
//------------------------------------------------------------------------------
}
#endif // PD_ELEMENT_H
