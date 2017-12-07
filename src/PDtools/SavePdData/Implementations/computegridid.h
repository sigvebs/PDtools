#ifndef COMPUTEGRIDID_H
#define COMPUTEGRIDID_H

#include "PDtools/SavePdData/savepddata.h"

//------------------------------------------------------------------------------
namespace PDtools {
class Grid;

class ComputeGridId : public ComputeProperty {
public:
  ComputeGridId(PD_Particles &particles, Grid &grid);
  ~ComputeGridId();
  virtual void update(const int id_i, const int i);

protected:
  Grid &m_grid;
  arma::mat &m_r;
  arma::mat &m_data;
  int m_indexgrid;
};
//------------------------------------------------------------------------------
}

#endif // COMPUTEGRIDID_H
