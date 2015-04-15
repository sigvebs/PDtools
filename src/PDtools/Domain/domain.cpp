#include "domain.h"

//------------------------------------------------------------------------------
PDtools::Domain::Domain(int dim, vector <pair<double, double>> boundaries):
    dim(dim), m_boundaries(boundaries)
{
    double x_len = m_boundaries[0].second - m_boundaries[0].first;
    double y_len = m_boundaries[1].second - m_boundaries[1].first;
    double z_len = m_boundaries[2].second - m_boundaries[2].first;

    m_boundaryLength = {x_len, y_len, z_len};
    m_periodicBoundaries = {0, 0, 0};
}
//------------------------------------------------------------------------------
void PDtools::Domain::setBoundaries(vector <pair<double, double>> boundaries)
{
    m_boundaries = boundaries;

    double x_len = m_boundaries[0].second - m_boundaries[0].first;
    double y_len = m_boundaries[1].second - m_boundaries[1].first;
    double z_len = m_boundaries[2].second - m_boundaries[2].first;

    m_boundaryLength[0] = x_len;
    m_boundaryLength[1] = y_len;
    m_boundaryLength[2] = z_len;
}
//------------------------------------------------------------------------------
