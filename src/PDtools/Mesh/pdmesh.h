#ifndef PDMESH_H
#define PDMESH_H

#include "config.h"

namespace PDtools {
//------------------------------------------------------------------------------
class PdMesh {
public:
  PdMesh();
  PdMesh(vector<array<double, 3>> vertices, vector<int> nodeIds,
         vector<vector<int>> elements, vector<int> elementIds);

  size_t nVertices() const;
  size_t nElements() const;

  const vector<array<double, 3>> &vertices() const;
  const vector<vector<int>> &elements() const;

  size_t nTriangles() const;
  size_t nQuads() const;

  unordered_map<int, int> getIdToCol_vertices() const;

  unordered_map<int, int> getIdToCol_elements() const;

private:
  size_t m_nVertices;
  size_t m_nElements;
  vector<array<double, 3>> m_vertices;
  vector<vector<int>> m_elements;
  vector<int> m_nodeIds;
  vector<int> m_elementIds;
  unordered_map<int, int> m_idToCol_nodes;
  unordered_map<int, int> m_idToCol_elements;

  size_t m_nTriangles;
  size_t m_nQuads;
};

//------------------------------------------------------------------------------
}
#endif // PDMESH_H
