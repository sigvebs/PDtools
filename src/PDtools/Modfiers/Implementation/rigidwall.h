#ifndef RIGIDWALL_H
#define RIGIDWALL_H

#include <unordered_map>
#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
class Grid;
class GridPoint;
//------------------------------------------------------------------------------
class RigidWall: public Modifier
{
public:
    RigidWall(Grid &grid, int orientationAxis, int topOrBottom, double lc);
    virtual void initialize();
    virtual void evaluateStepOne();

private:
    Grid &m_grid;
    int m_orientationAxis;
    int m_plane;
    vector<int> otherAxis;
    vector<GridPoint*> m_gridpoints;
    vector<pair<double, double> > m_boundary;
    int m_top = 0;
    int m_bottom = 0;
    double m_lc;
};
//------------------------------------------------------------------------------
}
#endif // RIGIDWALL_H
