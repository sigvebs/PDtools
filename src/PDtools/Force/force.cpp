#include "force.h"

#include <stdio.h>
#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
Force::Force(PD_Particles &particles):
    m_particles(particles),
    m_r(m_particles.r()),
    m_r0(m_particles.r0()),
    m_F(m_particles.F()),
    m_data(m_particles.data()),
    m_pIds(m_particles.pIds())
{

}
//------------------------------------------------------------------------------
Force::~Force()
{

}
//------------------------------------------------------------------------------
void Force::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    m_E = E;
    m_nu = nu;
    m_h = h;
    m_delta = delta;
    m_dim = dim;
    m_lc = lc;
}
//------------------------------------------------------------------------------
double Force::calculatePotentialEnergyDensity(const std::pair<int, int> &idCol)
{
    (void) idCol;
    return 0;
}
//------------------------------------------------------------------------------
void Force::calculatePotentialEnergy(const std::pair<int, int> &idCol,
                                     int indexPotential)
{
    (void) idCol;
    (void) indexPotential;
    //    std::cerr << "ERROR: potential energy not implemented for this force" << std::endl;
}
//------------------------------------------------------------------------------
double Force::calculateBondEnergy(const std::pair<int, int> &idCol,
                                  std::pair<int, std::vector<double> > &con)
{
    (void) idCol;
    (void) con;

    return 0;
}
//------------------------------------------------------------------------------
void Force::calculateStress(const std::pair<int, int> &idCol, const int (&indexStress)[6])
{
    (void) idCol;
    (void) indexStress;
    //    std::cerr << "ERROR: stress not implemented for this force" << std::endl;
}
//------------------------------------------------------------------------------
void Force::updateState()
{

}
//------------------------------------------------------------------------------
void Force::updateState(const std::pair<int, int> &idCol)
{
    (void) idCol;
}
//------------------------------------------------------------------------------
double Force::calculateStableMass(const std::pair<int, int> &idCol, double dt)
{
    (void) idCol;
    (void) dt;

    return 0.;
}
//------------------------------------------------------------------------------
void Force::numericalInitialization(bool ni)
{
    m_numericalInitialization = ni;
}
//------------------------------------------------------------------------------
void Force::setDim(int dim)
{
    m_dim = dim;
}
//------------------------------------------------------------------------------
void Force::applySurfaceCorrection(double strain)
{
    applyStrainCorrection(strain);
//    applyShearCorrection(strain);
}
//------------------------------------------------------------------------------
void Force::applyStrainCorrection(double strain)
{
    arma::vec3 scaleFactor;
    arma::mat g = arma::zeros(m_particles.nParticles(), m_dim);

    const int iDr0 = m_particles.getPdParamId("dr0");
    const int iForceScaling = m_particles.getPdParamId("forceScalingBond");

    // Stretching all particle in the x-direction
    scaleFactor(0) = strain;
    scaleFactor(1) = 0;
    scaleFactor(2) = 0;
//    scaleFactor(1) = -m_nu*strain;
//    scaleFactor(2) = -m_nu*strain;

    double W_infty = 0;
    switch(m_dim)
    {
    case 3:
        W_infty = 0.6*m_E*pow(strain, 2);
        break;
    case 2:
        W_infty = 9./16.*m_E*pow(strain, 2);
        break;
    case 1:
        W_infty = 0.5*m_E*pow(strain, 2);
        break;
    }

    for(int a=0; a<m_dim; a++)
    {
        if(a == 1)
            scaleFactor.swap_rows(0,1);
        else if(a == 2)
            scaleFactor.swap_rows(1,2);

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Loading the geometry
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                m_r(col_i, d) = (1 + scaleFactor(d))*m_r(col_i, d);
            }
        }

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            double W = this->calculatePotentialEnergyDensity(idCol);
            g(col_i, a) = W;
        }

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
        // Resetting the positions
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                m_r(col_i, d) = m_r(col_i, d)/(1 + scaleFactor(d));
            }
        }
    }

    // Scaling the energy with the median energy, which we assume
    // to be the bulk energy
    for(int a=0; a<m_dim; a++)
    {
        //#ifdef USE_OPENMP
        //# pragma omp parallel for
        //#endif
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            const double g_i = g(col_i, a);
            const double W =  W_infty/g_i;
            g(col_i, a) = W;
        }
    }

    // Calculating the scaling
//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        const int pId = idCol.first;
        const int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = m_pIds[id_j];

            const double dr0Len = con.second[iDr0];
            const arma::vec3 &n = (m_r.row(col_i).t() - m_r.row(col_j).t())/dr0Len;

            arma::vec3 g_mean;
            double G = 0;
            for(int d=0; d<m_dim; d++)
            {
                g_mean(d) = 0.5*(g(col_i, d) + g(col_j, d));
                G += pow(n(d)/g_mean(d), 2);
            }

            G = pow(G, -0.5);
            con.second[iForceScaling] *= G;
        }
    }
}
//------------------------------------------------------------------------------
void Force::applyShearCorrection(double shear)
{
    double m_mu = 0.5*m_E/(1 + m_nu); // Tmp

    arma::vec3 strainFactor;
    arma::vec3 scaleFactor;
    arma::mat g = arma::zeros(m_particles.nParticles(), m_dim);

    const int iDr0 = m_particles.getPdParamId("dr0");
    const int iForceScaling = m_particles.getPdParamId("forceScalingBond");

    // Performing a simple shear of all particle in the x, y and z-direction
    arma::ivec3 axis;
    strainFactor(0) = shear;
    strainFactor(1) = 0;
    strainFactor(2) = 0;
    axis(0) = 1;
    axis(1) = 0;
    axis(2) = 0;

    for(int a=0; a<m_dim; a++)
    {
        if(a == 1)
        {
            strainFactor.swap_rows(1,2);
            strainFactor.swap_rows(0,1);
            axis(0) = 2;
            axis(1) = 0;
            axis(2) = 1;
        }
        else if(a == 2)
        {
            strainFactor.swap_rows(2,0);
            strainFactor.swap_rows(1,2);
            axis(0) = 2;
            axis(1) = 0;
            axis(2) = 0;
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Loading the geometry
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                double shearFactor = strainFactor(d)*m_r(col_i, axis(d));
                m_r(col_i, d) = m_r(col_i, d) + shearFactor;
            }
        }

        double W_s = 0.5*m_mu*shear*shear;
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Calculating the elastic energy density
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            const int col_i = idCol.second;
            double W = this->calculatePotentialEnergyDensity(idCol);
            const double factor =  W_s/W;
            g(col_i, a) = factor;
        }

#ifdef USE_OPENMP
# pragma omp parallel for
#endif
        // Resetting the positions
        for(unsigned int i=0; i<m_particles.nParticles(); i++)
        {
            pair<int, int> idCol(i, i);
            int col_i = idCol.second;

            for(int d=0; d<m_dim; d++)
            {
                m_r(col_i, d) = m_r0(col_i, d); // Tmp
            }
        }
    }

    // Calculating the scaling
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        pair<int, int> idCol(i, i);
        const int pId = idCol.first;
        const int col_i = idCol.second;

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);

        for(auto &con:PDconnections)
        {
            const int id_j = con.first;
            const int col_j = m_pIds[id_j];

            const double dr0Len = con.second[iDr0];
            const arma::vec3 &n = (m_r.row(col_i).t() - m_r.row(col_j).t())/dr0Len;

            arma::vec3 g_mean;
            double G = 0;
            for(int d=0; d<m_dim; d++)
            {
                g_mean(d) = 0.5*(g(col_i, d) + g(col_j, d));
                G += pow(n(d)/g_mean(d), 2);
            }

            G = pow(G, -0.5);
            con.second[iForceScaling] *= G;
        }
    }
}
//------------------------------------------------------------------------------
}
