#ifndef PDFUNCTIONS_H
#define PDFUNCTIONS_H

#include <array>
#include <vector>
#include <armadillo>
#include "pdfunctionsmpi.h"

using namespace std;
//------------------------------------------------------------------------------
namespace PDtools
{
class PD_Particles;
class Grid;
class Force;

void
setPdConnections(PD_Particles & particles, Grid & grid,
                      double radius, double lc);

void
addFractures(PD_Particles & particles, const vector<pair<double,double>> &domain, const vector<double> &fracture);

void
removeVoidConnections(PD_Particles & particles, Grid & grid,
                      double delta, double lc);

void
cleanUpPdConnections(PD_Particles & particles);

void
applyVolumeCorrection(PD_Particles &particles, double delta, double lc, int dim);

void
reCalculatePdMicromodulus(PD_Particles & particles, int dim);

void
reCalculatePdFractureCriterion(PD_Particles & particles, double G0,
                                    double delta, double h=-1);
void
calculateRadius(PD_Particles & particles, int dim, double h=1.);

void
surfaceCorrection(PD_Particles & particles, vector<Force*> &forces, double E, double nu, int dim);

void
applyInitialStrainStrain(PD_Particles & particles, double strain, int axis, pair<double, double>area);

void
setPD_N3L(PD_Particles & particles);

vector<array<int, 3>>
possibleConfigurations2d(const int n);

vector<array<int, 3>>
possibleConfigurations3d(const int n);

vector<int>
optimalConfigurationCores(const int nCores, const vector<double> &domain, const int dim);

string
getFileEnding(string filename);
}
//------------------------------------------------------------------------------
#endif // PDFUNCTIONS_H
