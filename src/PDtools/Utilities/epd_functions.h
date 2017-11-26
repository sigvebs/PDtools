#ifndef EPD_FUNCTIONS_H
#define EPD_FUNCTIONS_H


#include <fstream>
using namespace std;

namespace PDtools
{
class PD_Particles;
class Grid;

//------------------------------------------------------------------------------
int saveElementQuadrature(PD_Particles &nodes, PDtools::Grid &grid, const string &savePath);
//------------------------------------------------------------------------------
}
#endif // EPD_FUNCTIONS_H
