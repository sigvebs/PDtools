#include "pdmesh.h"
//------------------------------------------------------------------------------
namespace PDtools {
PdMesh::PdMesh() {}
//------------------------------------------------------------------------------
PdMesh::PdMesh(vector<array<double, 3>> _vertices, vector<int> nodeIds,
               vector<vector<int>> _elements, vector<int> elementIds)
    : m_vertices(_vertices), m_elements(_elements), m_nodeIds(nodeIds),
      m_elementIds(elementIds) {
  m_nVertices = m_vertices.size();
  m_nElements = m_elements.size();

  for (size_t i = 0; i < m_nVertices; i++) {
    const int id_i = m_nodeIds[i];
    m_idToCol_nodes[id_i] = i;
  }

  for (size_t i = 0; i < m_nElements; i++) {
    const int id_i = m_elementIds[i];
    m_idToCol_elements[id_i] = i;
  }

  m_nTriangles = 0;
  m_nQuads = 0;

  for (const vector<int> &element : m_elements) {
    switch (element.size()) {
    case 3:
      m_nTriangles++;
      break;
    case 4:
      m_nQuads++;
      break;
    default:
      cerr << "Mesh containts polygons of non triangle/quad shape" << endl;
      exit(EXIT_FAILURE);
    }
  }
}
//------------------------------------------------------------------------------
size_t PdMesh::nVertices() const { return m_nVertices; }
//------------------------------------------------------------------------------
size_t PdMesh::nElements() const { return m_nElements; }
//------------------------------------------------------------------------------
const vector<array<double, 3>> &PdMesh::vertices() const { return m_vertices; }
//------------------------------------------------------------------------------
const vector<vector<int>> &PdMesh::elements() const { return m_elements; }
//------------------------------------------------------------------------------
size_t PdMesh::nTriangles() const { return m_nTriangles; }
//------------------------------------------------------------------------------
size_t PdMesh::nQuads() const { return m_nQuads; }
//------------------------------------------------------------------------------
unordered_map<int, int> PdMesh::getIdToCol_vertices() const {
  return m_idToCol_nodes;
}
//------------------------------------------------------------------------------
unordered_map<int, int> PdMesh::getIdToCol_elements() const {
  return m_idToCol_elements;
}
//------------------------------------------------------------------------------
}
