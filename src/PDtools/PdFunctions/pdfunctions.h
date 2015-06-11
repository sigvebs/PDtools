#ifndef PDFUNCTIONS_H
#define PDFUNCTIONS_H

#include <vector>

using namespace std;
//------------------------------------------------------------------------------
namespace PDtools
{
class PD_Particles;
class Grid;
class Force;

void setPdConnections(PD_Particles & particles, Grid & grid,
                      double radius, double alpha=0.25);
void reCalculatePdMicromodulus(PD_Particles & particles, int dim);
void reCalculatePdFractureCriterion(PD_Particles & particles, double G0,
                                    double delta, double h=-1);
void calculateRadius(PD_Particles & particles, int dim, double h=1.);
void surfaceCorrection(PD_Particles & particles, vector<Force*> &forces, double k, double nu, int dim);
}
//------------------------------------------------------------------------------
#endif // PDFUNCTIONS_H
