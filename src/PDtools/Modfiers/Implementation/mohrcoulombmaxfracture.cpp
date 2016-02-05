#include "mohrcoulombmaxfracture.h"

#include "PDtools/Particles/pd_particles.h"

#include <stdlib.h>

namespace PDtools
{
//------------------------------------------------------------------------------
MohrCoulombMaxFracture::MohrCoulombMaxFracture(double mu, double C, double T, int dim):
    m_C(C), m_T(T), m_dim(dim)
{
    m_phi = mu*M_PI/180.;
    m_d = tan(m_phi);
    m_neededProperties= {pair<string, int>("stress", 1)};
    m_C = 0.5*C*(1./(sqrt(m_d*m_d + 1.) + m_d));
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::registerParticleParameters()
{
    m_data = &m_particles->data();
    m_indexUnbreakable =  m_particles->registerParameter("unbreakable");
    m_indexConnected = m_particles->registerPdParameter("connected");
    m_indexCompute = m_particles->registerPdParameter("compute");
    m_indexBroken = m_particles->registerParameter("broken", 0);
    m_idToCol = &m_particles->idToCol();

    switch(m_dim)
    {
    case 1:
        m_ghostParameters = {"nx"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexNormal[0] = m_particles->registerParameter("nx");
        break;
    case 2:
        m_ghostParameters = {"nx", "ny"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");

        m_indexNormal[0] = m_particles->registerParameter("nx");
        m_indexNormal[1] = m_particles->registerParameter("ny");
        break;
    case 3:
        m_ghostParameters = {"nx", "ny", "nz"};
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStress[3] = m_particles->registerParameter("s_zz");
        m_indexStress[4] = m_particles->registerParameter("s_xz");
        m_indexStress[5] = m_particles->registerParameter("s_yz");
        m_indexNormal[0] = m_particles->registerParameter("nx");
        m_indexNormal[1] = m_particles->registerParameter("ny");
        m_indexNormal[2] = m_particles->registerParameter("nz");
        break;
    }
    m_ghostParameters.push_back("broken");
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::initialize()
{
    srand (time(NULL));
    m_broken = false;
    m_cosTheta = cos(M_PI/2. + m_phi);
    m_sinTheta = sin(M_PI/2. + m_phi);
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepOne(const int id_i, const int i)
{
    mat & data = *m_data;
    const mat & r = m_particles->r();

    if(data(i, m_indexUnbreakable) >= 1)
        return;

    const int broken_i = data(i, m_indexBroken);
    double n_i[m_dim];
    double n_j[m_dim];
    double dr_ij[m_dim];

    if(broken_i)
    {
        for(int d=0; d<m_dim; d++)
        {
            n_i[d] = data(i, m_indexNormal[d]);
        }
    }

    vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

    for(auto &con:PDconnections)
    {
        const int id_j = con.first;
        const int j = (*m_idToCol)[id_j];

        if(data(j, m_indexUnbreakable) >= 1)
            continue;

        if(con.second[m_indexConnected] <= 0.5)
            continue;

        if(broken_i)
        {
            double dotproduct = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(j, d) - r(i, d);
                dotproduct += n_i[d]*dr_ij[d];
            }

            if(dotproduct > 0)
                con.second[m_indexConnected] = 0;
        }

        const int broken_j = data(j, m_indexBroken);
        if(broken_j)
        {
            double dotproduct = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = r(i, d) - r(j, d);
                n_j[d] = data(j, m_indexNormal[d]);
                dotproduct += n_j[d]*dr_ij[d];
            }

            if(dotproduct > 0)
                con.second[m_indexConnected] = 0;
        }
    }

}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepTwo(const int id_i, const int i)
{
    (void) id_i;
    mat & data = *m_data;

    if(m_dim == 2)
    {
        double sx, sy, sxy;
        sx = data(i, m_indexStress[0]);
        sy = data(i, m_indexStress[1]);
        sxy = data(i, m_indexStress[2]);

        double first = 0.5*(sx + sy);
        double second = sqrt(0.25*(sx - sy)*(sx - sy) + sxy*sxy);

        double p_2 = first + second; // max
        double p_1 = first - second; // min

        double shear = fabs(0.5*(p_1 - p_2)*m_sinTheta);
        double normal = 0.5*(p_1 + p_2) + 0.5*(p_1 - p_2)*m_cosTheta;

        double s_max = max(sx, sy);
        double s_min = min(sx, sy);
        double theta = 0.5*atan(2.*fabs(sxy)/(s_max - s_min));
        int broken = 0;

        if(shear >= fabs(m_C - m_d*normal) && normal < 0)
        {
            data(i, m_indexBroken) = 1;
            m_broken = true;
            theta += M_PI/4.;
            broken = 1;
            cout << id_i << " " << theta*180./M_PI << " n = ";
            cout << data(i, m_indexNormal[0]) << " " << data(i, m_indexNormal[1]) << endl;
        }
        else if(p_2 >= m_T)
        {
            data(i, m_indexBroken) = 1;
            m_broken = true;
            broken = 1;
        }
        else
        {
            data(i, m_indexBroken) = 0;
        }

        if(broken)
        {
            int ra = rand() % 2;
            if(ra >= 1)
                theta += M_PI;

            if(sy < sx)
            {
                data(i, m_indexNormal[0]) = cos(theta);
                data(i, m_indexNormal[1]) = sin(theta);
            }
            else
            {
                data(i, m_indexNormal[0]) = -sin(theta);
                data(i, m_indexNormal[1]) = cos(theta);
            }
        }
    }
}
//------------------------------------------------------------------------------
void MohrCoulombMaxFracture::evaluateStepTwo()
{
    if(m_broken)
    {
        m_state = true;
    }
    else
    {
        m_state = false;
    }

    m_broken = false;
}
//------------------------------------------------------------------------------
}
