#include "calculatestrain.h"

#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculateStrain::CalculateStrain(double delta, vector<pair<double,double>> &domain):
    CalculateProperty("strain"),
    m_delta(delta),
    m_domain(domain)
{
}
//------------------------------------------------------------------------------
void CalculateStrain::initialize()
{
    m_indexConnected = m_particles->getPdParamId("connected");
    m_indexVolume = m_particles->getParamId("volume");
    m_indexVolumeScaling = m_particles->getPdParamId("volumeScaling");
    m_indexDr0 = m_particles->getPdParamId("dr0");

    if(m_dim == 2)
        m_h = m_domain[2].second - m_domain[2].first;

    switch(m_dim)
    {
    case 1:
        m_nStrainElements = 1;
        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        break;
    case 2:
        m_nStrainElements = 3;
        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        m_indexStrain[1] = m_particles->registerParameter("e_yy");
        m_indexStrain[2] = m_particles->registerParameter("e_xy");
        break;
    case 3:
        m_nStrainElements = 6;
        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        m_indexStrain[1] = m_particles->registerParameter("e_yy");
        m_indexStrain[2] = m_particles->registerParameter("e_xy");
        m_indexStrain[3] = m_particles->registerParameter("e_zz");
        m_indexStrain[4] = m_particles->registerParameter("e_xz");
        m_indexStrain[5] = m_particles->registerParameter("e_yz");
        break;
    }
}
//------------------------------------------------------------------------------
void CalculateStrain::clean()
{
    const int nParticles = m_particles->nParticles();
    arma::mat & data = m_particles->data();

    // Zeroing the stress
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++)
    {
        for(int s=0; s<m_nStrainElements; s++)
        {
            data(i, m_indexStrain[s]) = 0;
        }
    }
}
//------------------------------------------------------------------------------
void CalculateStrain::update()
{
    const ivec &colToId = m_particles->colToId();
    const auto &idToCol = m_particles->idToCol();
    const int nParticles = m_particles->nParticles();
    const mat &r = m_particles->r();
    const mat &r0 = m_particles->r0();
    const mat &data = m_particles->data();

    mat K = zeros(m_dim, m_dim);
    double dr_ij[m_dim];
    double dr0_ij[m_dim];

    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId(i);
        vector<pair<int, vector<double>>> & PDconnections_i = m_particles->pdConnections(id_i);
        const int nConnections = PDconnections_i.size();
        K.zeros();

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            const auto &con_i = PDconnections_i[l_j];
            if(con_i.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con_i.first;
            const int j = idToCol.at(id_j);

            const double vol_j = data(j, m_indexVolume);
            const double dr0 = con_i.second[m_indexDr0];
            const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
            const double v = vol_j*volumeScaling_ij;
            const double w = m_delta/dr0;

            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(j, d) - r(i, d);
                dr0_ij[d] = r0(j, d) - r0(i, d);
            }

//            kxx += w*dr0_ij[0]*dr0_ij[0]*v;

            for(int d=0; d<m_dim; d++)
            {
                for(int d2=0; d2<m_dim; d2++)
                {
                    K(d, d2) += w*dr0_ij[d]*dr0_ij[d2]*v;
                }
            }
        }
//        double fac = (pow(m_delta, 4)*M_PI*m_h)/3.;

    }
}
//------------------------------------------------------------------------------
}
