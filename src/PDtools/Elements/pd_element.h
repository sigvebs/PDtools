#ifndef PD_ELEMENT_H
#define PD_ELEMENT_H

#include "config.h"

//------------------------------------------------------------------------------
namespace PDtools {
class PdMesh;
class PD_Particles;
class Grid;

//------------------------------------------------------------------------------
template <size_t T> class PD_element {
public:
  PD_element() {}
  PD_element(size_t _id, array<size_t, T> _verticeIds)
      : m_id(_id), m_verticeIds(_verticeIds) {}
  const array<size_t, T> &verticeIds() const { return m_verticeIds; }

  size_t id() const { return m_id; }

  const arma::mat &guassianQuadraturePoints_initial() const {
    return m_guassianQuadraturePoints_initial;
  }
  void setGuassianQuadraturePoints_initial(
      const arma::mat &guassianQuadraturePoints_initial) {
    m_guassianQuadraturePoints_initial =
        arma::mat(guassianQuadraturePoints_initial);
  }

  arma::mat &guassianQuadraturePoints() { return m_guassianQuadraturePoints; }
  void setGuassianQuadraturePoints(const arma::mat &guassianQuadraturePoints) {
    m_guassianQuadraturePoints = arma::mat(guassianQuadraturePoints);
  }

  arma::mat guassianQuadratureWeights() const {
    return m_guassianQuadratureWeights;
  }
  void
  setGuassianQuadratureWeights(const arma::vec &guassianQuadratureWeights) {
    m_guassianQuadratureWeights = arma::vec(guassianQuadratureWeights);
  }

protected:
  size_t m_id;
  array<size_t, T> m_verticeIds;
  arma::mat m_guassianQuadraturePoints_initial;
  arma::mat m_guassianQuadraturePoints;
  arma::vec m_guassianQuadratureWeights;
};

class PD_lineElement : public PD_element<2> {
public:
  PD_lineElement() {}
  PD_lineElement(size_t _id, array<size_t, 2> _verticeIds)
      : PD_element(_id, _verticeIds) {}
};

class PD_triElement : public virtual PD_element<3> {
public:
  PD_triElement() {}
  PD_triElement(size_t _id, array<size_t, 3> _verticeIds)
      : PD_element(_id, _verticeIds) {}
};

class PD_quadElement : public PD_element<4> {
public:
  PD_quadElement() {}
  PD_quadElement(size_t _id, array<size_t, 4> _verticeIds)
      : PD_element(_id, _verticeIds) {}
};

// typedef PD_element<2> PD_lineElement;
// typedef PD_element<3> PD_triElement;
// typedef PD_element<4> PD_quadElement;

PD_Particles initializeElementPd(const PdMesh &msh, const Grid &grid,
                                 const size_t quadratureDegree);
void updateElementQuadrature(PD_Particles &discretization);
void printElement(const int element_id, PD_Particles &nodes);
//------------------------------------------------------------------------------
// class PD_elementsAndNodes
//{
// public:
//    PD_elementsAndNodes(Grid & grid);

//    void initialize(const PdMesh &msh);

// protected:
//    const int m_dim;
//    const Grid &m_grid;
//    PD_Particles m_particles;
//    vector<PD_triElement> m_triElements;
//    vector<PD_quadElement> m_quadElements;
//};
//------------------------------------------------------------------------------
}
#endif // PD_ELEMENT_H
