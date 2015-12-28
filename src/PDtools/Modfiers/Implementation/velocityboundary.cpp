#include "velocityboundary.h"

#include "PDtools/Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
VelocityBoundary::VelocityBoundary(double velAmplitude,
                                   double velOrientation,
                                   pair<double, double> boundary,
                                   int boundaryOrientation,
                                   double dt,
                                   int steps, bool isStatic)
{
    m_velAmplitude = velAmplitude;
    m_velOritentation = velOrientation;
    m_boundary = boundary;
    m_boundaryOrientation = boundaryOrientation;
    m_dv = velAmplitude / steps;
    m_v = 0;
    m_dt = dt;
    m_isStatic = isStatic;

    if(m_velOritentation == 0)
    {
        m_otherAxis = {1, 2};
    }
    else if(m_velOritentation == 1)
    {
        m_otherAxis = {0, 2};
    }
    else if(m_velOritentation == 2)
    {
        m_otherAxis = {1, 2};
    }
}
//------------------------------------------------------------------------------
VelocityBoundary::~VelocityBoundary()
{

}
//------------------------------------------------------------------------------
void VelocityBoundary::evaluateStepOne()
{
    if(fabs(m_v) < fabs(m_velAmplitude))
    {
        m_v += m_dv;
    }

    arma::mat & v = m_particles->v();
    arma::mat & r = m_particles->r();
    arma::mat & F = m_particles->F();
//    arma::mat &data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();

//    const int colRho = m_particles->getParamId("rho");
    const double v_dt = m_v * m_dt;
//    const double dtRhoHalf = 0.5*m_dt;

    for(pair<int, int> &idCol:m_boundaryParticles)
    {
        int col_i = idCol.second;

        v(col_i, m_boundaryOrientation) = m_v;
        if(isStatic(col_i))
            r(col_i, m_boundaryOrientation) += v_dt;

//        double rho = data(col_i, colRho);
//        double dtRho = dtRhoHalf/rho;
//        for(int d:m_otherAxis)
//        {
//            v(col_i, d) += F(col_i, d)*dtRho;
//            r(col_i, d) += v(col_i, d)*m_dt;
//        }
//        for(int d:m_otherAxis)
//        {
//            F(col_i, d) = 0;
//        }
        F(col_i, m_boundaryOrientation) = 0;
    }
}
//------------------------------------------------------------------------------
void VelocityBoundary::evaluateStepTwo()
{
//    arma::mat & v = m_particles->v();
//    arma::mat & F = m_particles->F();
//    arma::mat &data = m_particles->data();
//    const int colRho = m_particles->getParamId("rho");
//    const double dtRhoHalf = 0.5*m_dt;


//    for(pair<int, int> &idCol:m_boundaryParticles)
//    {
//        int col_i = idCol.second;

//        double rho = data(col_i, colRho);
//        double dtRho = dtRhoHalf/rho;
//        for(int d:m_otherAxis)
//        {
//            v(col_i, d) += F(col_i, d)*dtRho;
//        }
//    }
}
//------------------------------------------------------------------------------
void VelocityBoundary::initialize()
{
    // Selecting particles
    const arma::mat & r = m_particles->r();
    arma::mat & data = m_particles->data();
    arma::imat & isStatic = m_particles->isStatic();
    int unbreakablePos;

    if(m_particles->hasParameter("unbreakable"))
    {
        unbreakablePos = m_particles->getParamId("unbreakable");
    }
    else
    {
        unbreakablePos = m_particles->registerParameter("unbreakable");
    }

    double uRadius = 0.5*(m_boundary.second - m_boundary.first);
    int n = 0;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        int col_i = i;
        double pos = r(col_i, m_boundaryOrientation);
        if(m_boundary.first <= pos && pos < m_boundary.second)
        {
            pair<int, int> pId(i, i);
            if(m_isStatic)
                isStatic(pId.second) = 1;
#ifdef USE_OPENMP
#pragma omp critical
#endif
            m_boundaryParticles.push_back(pId);
            data(i, unbreakablePos) = 1;
        }

        if(m_boundary.first - uRadius <= pos && pos < uRadius + m_boundary.second)
        {
//#pragma omp critical
//            data(i, unbreakablePos) = 1;
            n++;
        }
    }
}
//------------------------------------------------------------------------------
}
