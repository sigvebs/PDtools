#include "epd_bondforce.h"
#include "PDtools/Particles/pd_particles.h"

#define OVERLAP 0
#define USE_VISCOUS 0
#define USE_CORRECTION 1

namespace PDtools
{
//------------------------------------------------------------------------------
EPD_bondForce::EPD_bondForce(PD_Particles &particles):
    Force(particles, "EPD bond force")
{
    m_indexMicromodulus = m_particles.registerParameter("micromodulus", 1);
    m_iOverlap = m_particles.getPdParamId("overlap");
    m_iConnected = m_particles.getPdParamId("connected");

#if USE_VISCOUS
    particles.setNeedGhostVelocity(true);
#endif

    //    m_ghostParameters = {"volume", "micromodulus"};
    //    m_initialGhostParameters = {"volume", "micromodulus"};
}
//------------------------------------------------------------------------------
void EPD_bondForce::calculateForces(const int id_i, const int i)
{
    const double c = m_data(i, m_indexMicromodulus);
    vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

    const int nConnections = PDconnections.size();
    double dr_ij[m_dim];
    double dr0_ij[m_dim];

#if USE_CORRECTION
    vec f_correction = zeros(m_dim);
    vec du_ij = zeros(m_dim);
#endif
    for(int l_j=0; l_j<nConnections; l_j++)
    {
        auto &con = PDconnections[l_j];
        if(con.second[m_iConnected] <= 0.5)
            continue;

        const int polygon_id = con.first;
        const double overlap = con.second[m_iOverlap];
        const int polygon_i = m_idToElement.at(polygon_id);

        PD_quadElement & element = m_quadElements[polygon_i];
        const mat & gaussPoints_initial = element.guassianQuadraturePoints_initial();
        const mat & gaussPoints = element.guassianQuadraturePoints();
        const vec & gaussWeights = element.guassianQuadratureWeights();

        const int nIntegrationPoints = gaussPoints.n_rows;

        for(int j=0; j<nIntegrationPoints; j++) {
            const double w_j = gaussWeights[j];

            double dr2 = 0;
            double dr2_0 = 0;
            for(int d=0; d<m_dim; d++)
            {
                dr_ij[d] = gaussPoints(j, d) - m_r(i, d);
                dr2 += dr_ij[d]*dr_ij[d];
                dr0_ij[d] = gaussPoints_initial(j, d) - m_r0(i, d);
                dr2_0 += dr0_ij[d]*dr0_ij[d];
#if USE_CORRECTION
                du_ij[d] =  dr_ij[d] - dr0_ij[d];
#endif
            }

#if OVERLAP
            const double dr0 = sqrt(dr2_0);
            const double dr = sqrt(dr2);
            const double ds = dr - dr0;
            if(ds < 1e-12)
                continue;
            const double s = ds/dr0;
            const double weights = weightFunction(dr0);
            const double fbond_ij = overlap*weights*c*s*w_j/dr;
#else
            if(dr2_0> m_delta*m_delta)
                continue;
            const double dr0 = sqrt(dr2_0);
            const double dr = sqrt(dr2);
            const double ds = dr - dr0;
            //            if(ds < 1e-12)
            //                continue;
            const double s = ds/dr0;
            const double weights = weightFunction(dr0);
            const double fbond_ij = weights*c*s*w_j/dr;

#if USE_CORRECTION
            f_correction += weights*m_dampCoeff*du_ij*w_j;
#endif

#endif
            for(int d=0; d<m_dim; d++)
            {
                const double Fd = dr_ij[d]*fbond_ij;
                m_F(i, d) += Fd;
            }
        }
    }
    for(int d=0; d<m_dim; d++) {
        m_F(i, d) += f_correction(d);
    }
}
//------------------------------------------------------------------------------
double EPD_bondForce::calculateStableMass(const int id_a, const int a, double dt)
{
    dt *= 1.1;

    const double c = m_data(a, m_indexMicromodulus);

    double m[m_dim];
    double dr0_ij[m_dim];

    for(int d=0; d<m_dim;d++)
    {
        m[d] = 0;
    }

    const vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_a);

    double k[m_dim];

    for(int i=0; i<m_dim; i++) {
        for(int d=0; d<m_dim;d++) {
            k[d] = 0;
        }

        for(auto &con:PDconnections) {
            if(con.second[m_iConnected] <= 0.5)
                continue;

            const int polygon_id = con.first;
            const int polygon_i = m_idToElement.at(polygon_id);

            PD_quadElement & element = m_quadElements[polygon_i];
            const mat & gaussPoints_initial = element.guassianQuadraturePoints_initial();
            const mat & gaussPoints = element.guassianQuadraturePoints();
            const vec & gaussWeights = element.guassianQuadratureWeights();

            const int nIntegrationPoints = gaussPoints.n_rows;
            for(int j=0; j<nIntegrationPoints; j++) {
                const double w_j = gaussWeights[j];

                double dr2_0 = 0;
                for(int d=0; d<m_dim; d++)
                {
                    dr0_ij[d] = gaussPoints_initial(j, d) - m_r0(a, d);
                    dr2_0 += dr0_ij[d]*dr0_ij[d];
                }

                if(dr2_0> m_delta*m_delta)
                    continue;
                const double dr0 = sqrt(dr2_0);

                const double weightFunc = weightFunction(dr0);
                const double dr0Len3 = pow(dr0, 3);
                double coeff = c*weightFunc*w_j/dr0Len3;
#if OVERLAP
                const double overlap = con.second[m_iOverlap];
                coeff *= overlap;
#endif

                double sum = 0;

                for(int d=0; d<m_dim; d++)
                {
                    sum += fabs(dr0_ij[d]);
                }
                sum *= fabs(dr0_ij[i])*coeff;

                k[i] += sum;
            }
        }
        m[i] = k[i];
    }

    double stiffness = 0;

    for(int d=0;d<m_dim; d++){
        if(m[d]>stiffness){
            stiffness = m[d];
        }
    }

    return 4.*0.25*pow(dt, 2)*stiffness;
}
//------------------------------------------------------------------------------
void EPD_bondForce::initialize(double E, double nu, double delta, int dim, double h, double lc)
{
    Force::initialize(E, nu, delta, dim, h, lc);
    double k;
    double c;

    m_dampCoeff = 3.*5.e-11; // 3*dt
    m_dampCoeff = 1.e13*E/(lc*lc);
    m_dampCoeff = 5.e13*E/(lc*lc); // Default value

    if(dim == 3)
    {
        nu = 1./4.;
        k = E/(3.*(1. - 2.*nu));
        c = 18.*k/(M_PI*pow(delta, 4));
    }
    else if(dim == 2)
    {
        nu = 1./3.;
        k = E/(2.*(1. - nu));
        c = 12*k/(h*M_PI*pow(delta, 3));
    }
    else if(dim == 1)
    {
        nu = 1./4.;
        k = E;
        c = 2*E/(h*h*pow(delta, 2));
    }
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported" << endl;
        cerr << "use 1, 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }

    if(m_numericalInitialization){
        m_particles.setParameter("micromodulus", k);
    }
    else
    {
        m_particles.setParameter("micromodulus", c);
    }

    computeMicromodulus();
}
//------------------------------------------------------------------------------
void EPD_bondForce::computeMicromodulus()
{
    // Scaling the bond force
    double nu = 1./3.;
    double k = m_E/(2.*(1. - nu));
    const double dimScaling = 2.*k*pow(m_dim, 2.);
    double dr0_ij[m_dim];
    const ivec & colToId = m_particles.colToId();
    double c_anal = 12*k/(m_h*M_PI*pow(m_delta, 3));

    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        const int id_i = colToId(i);
        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(id_i);

        const int nConnections = PDconnections.size();
        double m_i = 0;

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            auto &con = PDconnections[l_j];
            if(con.second[m_iConnected] <= 0.5)
                continue;

            const int polygon_id = con.first;
            const double overlap = con.second[m_iOverlap];
            const int polygon_i = m_idToElement.at(polygon_id);

            PD_quadElement & element = m_quadElements[polygon_i];
            const mat & gaussPoints_initial = element.guassianQuadraturePoints_initial();
            const vec & gaussWeights = element.guassianQuadratureWeights();

            const int nIntegrationPoints = gaussPoints_initial.n_rows;

            for(int j=0; j<nIntegrationPoints; j++) {
                const double w_j = gaussWeights[j];

                double dr2_0 = 0;
                for(int d=0; d<m_dim; d++)
                {
                    dr0_ij[d] = gaussPoints_initial(j, d) - m_r0(i, d);
                    dr2_0 += dr0_ij[d]*dr0_ij[d];
                }

#if OVERLAP
                const double dr0 = sqrt(dr2_0);
                m_i += overlap*dr0*weightFunction(dr0)*w_j;
#else
                const double dr0 = sqrt(dr2_0);
                if(dr0> m_delta)
                    continue;
                m_i += dr0*weightFunction(dr0)*w_j;
#endif
            }
        }

        const double c = m_data(i, m_indexMicromodulus);
        const double c_i = dimScaling/m_i;

        m_data(i, m_indexMicromodulus) =  0.9*c_i;
        //        m_data(i, m_indexMicromodulus) =  c;
    }
}
//------------------------------------------------------------------------------
}
