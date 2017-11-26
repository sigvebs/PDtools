#include "calculatestressstrain.h"

#include "Particles/pd_particles.h"
#include "PDtools/Force/force.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculateStressStrain::CalculateStressStrain(vector<Force *> &forces, double E,
                                             double nu, double delta,
                                             bool planeStress,
                                             bool greenLagrange):
    CalculateProperty("stress"),
    m_E(E),
    m_nu(nu),
    m_delta(delta),
    m_planeStress(planeStress)
{
    for(Force* force: forces) {
        m_forces.push_back(force);
    }
    m_greenStrain = greenLagrange;
}
//------------------------------------------------------------------------------
void CalculateStressStrain::initialize()
{
    m_particles->setNeedGhostR0(true);
    m_indexConnected = m_particles->getPdParamId("connected");
    m_indexVolume = m_particles->getParamId("volume");
    m_indexVolumeScaling = m_particles->getPdParamId("volumeScaling");
    m_indexDr0 = m_particles->getPdParamId("dr0");
    m_indexBrokenNow = m_particles->registerParameter("brokenNow", 0);

    switch(m_dim) {
    case 1:
        m_nStressStrainElements = 1;
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        m_indexK[0] = m_particles->registerParameter("K_00");
        break;
    case 2:
        m_nStressStrainElements = 3;
        m_indexK[0] = m_particles->registerParameter("K_00");
        m_indexK[1] = m_particles->registerParameter("K_11");
        m_indexK[2] = m_particles->registerParameter("K_01");
        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        m_indexStrain[1] = m_particles->registerParameter("e_yy");
        m_indexStrain[2] = m_particles->registerParameter("e_xy");
        break;
    case 3:
        m_nStressStrainElements = 6;
        m_indexK[0] = m_particles->registerParameter("K_00");
        m_indexK[1] = m_particles->registerParameter("K_11");
        m_indexK[2] = m_particles->registerParameter("K_01");
        m_indexK[3] = m_particles->registerParameter("K_22");
        m_indexK[4] = m_particles->registerParameter("K_02");
        m_indexK[5] = m_particles->registerParameter("K_12");

        m_indexStress[0] = m_particles->registerParameter("s_xx");
        m_indexStress[1] = m_particles->registerParameter("s_yy");
        m_indexStress[2] = m_particles->registerParameter("s_xy");
        m_indexStress[3] = m_particles->registerParameter("s_zz");
        m_indexStress[4] = m_particles->registerParameter("s_xz");
        m_indexStress[5] = m_particles->registerParameter("s_yz");

        m_indexStrain[0] = m_particles->registerParameter("e_xx");
        m_indexStrain[1] = m_particles->registerParameter("e_yy");
        m_indexStrain[2] = m_particles->registerParameter("e_xy");
        m_indexStrain[3] = m_particles->registerParameter("e_zz");
        m_indexStrain[4] = m_particles->registerParameter("e_xz");
        m_indexStrain[5] = m_particles->registerParameter("e_yz");
        break;
    }

    // Computing the shape tensor
    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();

    for(int i=0; i<nParticles; i++) {
        const int id_i = colToId(i);
        computeK(id_i, i);
    }

    m_lambda = m_E*m_nu/((1+m_nu)*(1. - 2.*m_nu));
    m_mu = 0.5*m_E/(1. + m_nu);
}
//------------------------------------------------------------------------------
void CalculateStressStrain::update()
{
    if(m_dim == 3) {
        update3d();
    } else if (m_dim == 2) {
        update2d();
    }
}
//------------------------------------------------------------------------------
void CalculateStressStrain::update2d()
{   const ivec &colToId = m_particles->colToId();
    const auto &idToCol = m_particles->idToCol();
    const int nParticles = m_particles->nParticles();
    const mat &r = m_particles->r();
    const mat &r0 = m_particles->r0();
    mat &data = m_particles->data();

    const double * d_volume =  data.colptr(m_indexVolume);
    double * d_indexBrokenNow =  data.colptr(m_indexBrokenNow);
    const double * m_x =  r.colptr(0);
    const double * m_y =  r.colptr(1);
    const double * m_z =  r.colptr(2);
    const double * m_x0 =  r0.colptr(0);
    const double * m_y0 =  r0.colptr(1);
    const double * m_z0 =  r0.colptr(2);


    mat F = zeros(m_dim, m_dim);
    mat K = zeros(m_dim, m_dim);
    mat P = zeros(m_dim, m_dim);
    mat strain = zeros(m_dim, m_dim);


    vector<double*> d_stress;
    vector<double*> d_strain;
    vector<double*> d_K;
    for(int i=0; i<m_nStressStrainElements; i++) {
        d_stress.push_back(data.colptr(m_indexStress[i]));
        d_strain.push_back(data.colptr(m_indexStrain[i]));
        d_K.push_back(data.colptr(m_indexK[i]));
    }

    double dr_ij[m_dim];
    double dr0_ij[m_dim];

    bool use_local_E = false;
    int iE = 0;
    if(m_particles->hasParameter("E_local")) {
        use_local_E = true;
        iE = m_particles->getParamId("E_local");
    }

    double E = m_E;
    double lambda = m_lambda;
    double mu = m_mu;
    const double nu = m_nu;

    const double E1 = 1./E;
    const double G = E/(2*(1+nu));
    const double G1 = 0.5/G;

    for(int i=0; i<nParticles; i++) {

#if 0
        P(0,0) = d_stress[0][i];
        P(1,1) = d_stress[1][i];
        P(0,1) = d_stress[2][i];

        if(m_planeStress) {
            d_strain[0][i] = E1*(P(0,0) - nu*P(1,1));
            d_strain[1][i] = E1*(P(1,1) - nu*P(0,0));
            d_strain[2][i] = E1*(1 + nu)*P(0,1);
        } else { // Plane strain
            d_strain[0][i] = E1*(1+nu)*((1-nu)*P(0,0) - nu*P(1,1));
            d_strain[1][i] = E1*(1+nu)*((1-nu)*P(1,1) - nu*P(0,0));
            d_strain[2][i] = E1*(1+nu)*P(0,1);
        }
#else
        const int id_i = colToId[i];

        if(d_indexBrokenNow[i]) {
            computeK(id_i, i);
            d_indexBrokenNow[i] = 0;
        }

        if(use_local_E) {
//            E = data(i, iE);
            lambda = E*nu/((1 + nu)*(1. - 2.*nu));
            mu = 0.5*E/(1. + nu);
        }

        K(0, 0) = d_K[0][i];
        K(1, 1) = d_K[1][i];
        K(0, 1) = d_K[2][i];
        K(1, 0) = K(0, 1);

        F.zeros();
        vector<pair<int, vector<double>>> & PDconnections_i = m_particles->pdConnections(id_i);
        const int nConnections = PDconnections_i.size();
        int nConnected = 0;

        for(int l_j=0; l_j<nConnections; l_j++) {

            const auto &con_i = PDconnections_i[l_j];
            if(con_i.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con_i.first;
            const int j = idToCol.at(id_j);

            const double vol_j = d_volume[j];
            const double dr0 = con_i.second[m_indexDr0];
            const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
            const double vol = vol_j*volumeScaling_ij;
            const double w = weightFunction(dr0);

            dr_ij[0] = m_x[j] - m_x[i];
            dr_ij[1] = m_y[j] - m_y[i];

            dr0_ij[0] = m_x0[j] - m_x0[i];
            dr0_ij[1] = m_y0[j] - m_y0[i];

            for(int d=0; d<m_dim; d++) {
                for(int d2=0; d2<m_dim; d2++) {
                    F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
                }
            }

            nConnected++;
        }

        if(nConnected <= 3) {
            for(int j=0; j<m_nStressStrainElements; j++) {
                data(i, m_indexStrain[j]) = 0;
                data(i, m_indexStress[j]) = 0;
            }
            continue;
        } else {
            F = F*K; // K = inv(K);
            if(m_greenStrain) {
                strain = 0.5*F.t()*F;
                for(int d=0; d<m_dim;d++) {
                    strain(d, d) -= 0.5;
                }
            } else {
                strain = 0.5*(F.t() + F);
                for(int d=0; d<m_dim;d++) {
                    strain(d, d) -= 1.;
                }
            }
        }

        data(i, m_indexStrain[0]) = strain(0,0);
        data(i, m_indexStrain[1]) = strain(1,1);
        data(i, m_indexStrain[2]) = strain(0,1);

        // Assuming linear elasticity
        // Constituent model, linear elastic
        // Computing the second PK stress
        if(m_planeStress) {
            const double a = E/(1. - nu*nu);
            P(0,0) = a*(strain(0,0) + nu*strain(1,1));
            P(1,1) = a*(strain(1,1) + nu*strain(0,0));
            P(0,1) = a*(1 - nu)*strain(0,1);
            P(1,0) = P(0,1);
        } else { // Plane strain
            const double a = E/((1. + nu)*(1. - 2*nu));
            P(0,0) = a*((1. - nu)*strain(0,0) + nu*strain(1,1));
            P(1,1) = a*((1. - nu)*strain(1,1) + nu*strain(0,0));
            P(0,1) = a*0.5*(1. - 2.*nu)*strain(0,1);
            P(1,0) = P(0,1);
        }

        // Converting PK2 stress to Cauchy stress
        if(m_greenStrain) {
            const double detF = 1./det(F);
            P = detF*F*P*F.t();
        }
        data(i, m_indexStress[0]) = P(0,0);
        data(i, m_indexStress[1]) = P(1,1);
        data(i, m_indexStress[2]) = P(0,1);
#endif
    }
}
//------------------------------------------------------------------------------
void CalculateStressStrain::update3d()
{
    const ivec &colToId = m_particles->colToId();
    const auto &idToCol = m_particles->idToCol();
    const int nParticles = m_particles->nParticles();
    const mat &r = m_particles->r();
    const mat &r0 = m_particles->r0();
    mat &data = m_particles->data();

    const double * d_volume =  data.colptr(m_indexVolume);
    double * d_indexBrokenNow =  data.colptr(m_indexBrokenNow);
    const double * m_x =  r.colptr(0);
    const double * m_y =  r.colptr(1);
    const double * m_z =  r.colptr(2);
    const double * m_x0 =  r0.colptr(0);
    const double * m_y0 =  r0.colptr(1);
    const double * m_z0 =  r0.colptr(2);

    mat F = zeros(m_dim, m_dim);
    mat K = zeros(m_dim, m_dim);
    mat P = zeros(m_dim, m_dim);
    mat strain = zeros(m_dim, m_dim);

    vector<double*> d_stress;
    vector<double*> d_strain;
    vector<double*> d_K;

    for(int i=0; i<m_nStressStrainElements; i++) {
        d_stress.push_back(data.colptr(m_indexStress[i]));
        d_strain.push_back(data.colptr(m_indexStrain[i]));
        d_K.push_back(data.colptr(m_indexK[i]));
    }

    double dr_ij[m_dim];
    double dr0_ij[m_dim];

    bool use_local_E = false;
    int iE = 0;
    int iNu = 0;
    if(m_particles->hasParameter("E_local")) {
        use_local_E = true;
        iE = m_particles->getParamId("E_local");
        iNu = m_particles->getParamId("nu_local");
    }

    double E = m_E;
    double lambda = m_lambda;
    double mu = m_mu;
    double nu = m_nu;

    double E1 = 1./E;
    double G = E/(2*(1+nu));
    double G1 = 0.5/G;

    for(int i=0; i<nParticles; i++) {
#if 0
        P(0,0) = d_stress[0][i];
        P(1,1) = d_stress[1][i];
        P(0,1) = d_stress[2][i];
        P(2,2) = d_stress[3][i];
        P(0,2) = d_stress[4][i];
        P(1,2) = d_stress[5][i];

        if(use_local_E) {
            E = data(i, iE);
            nu = data(i, iNu);
            E1 = 1./E;
            G = E/(2*(1+nu));
            G1 = 0.5/G;
        }

        d_strain[0][i] = E1*(P(0,0) - nu*(P(1,1) + P(2,2)));
        d_strain[1][i] = E1*(P(1,1) - nu*(P(0,0) + P(2,2)));
        d_strain[2][i] = G1*P(0,1);
        d_strain[3][i] = E1*(P(2,2) - nu*(P(1,1) + P(0,0)));
        d_strain[4][i] = G1*P(0,2);
        d_strain[5][i] = G1*P(1,2);
#else
        const int id_i = colToId[i];
        computeK(id_i, i);

        if(d_indexBrokenNow[i]) {
//            computeK(id_i, i);
            d_indexBrokenNow[i] = 0;
        }

        if(use_local_E) {
            //            E = data(i, iE);
            lambda = E*nu/((1 + nu)*(1. - 2.*nu));
            mu = 0.5*E/(1. + nu);
        }

        K(0, 0) = d_K[0][i];
        K(1, 1) = d_K[1][i];
        K(0, 1) = d_K[2][i];
        K(1, 0) = K(0, 1);
        K(2, 2) = d_K[3][i];
        K(0, 2) = d_K[4][i];
        K(1, 2) = d_K[5][i];
        K(2, 0) = K(0, 2);
        K(2, 1) = K(1, 2);

        F.zeros();
        vector<pair<int, vector<double>>> & PDconnections_i = m_particles->pdConnections(id_i);
        const int nConnections = PDconnections_i.size();
        int nConnected = 0;
        for(int l_j=0; l_j<nConnections; l_j++) {
            const auto &con_i = PDconnections_i[l_j];
            if(con_i.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con_i.first;
            const int j = idToCol.at(id_j);

            const double vol_j = d_volume[j];
            const double dr0 = con_i.second[m_indexDr0];
            const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
            const double vol = vol_j*volumeScaling_ij;
            const double w = weightFunction(dr0);

            dr_ij[0] = m_x[j] - m_x[i];
            dr_ij[1] = m_y[j] - m_y[i];
            dr_ij[2] = m_z[j] - m_z[i];

            dr0_ij[0] = m_x0[j] - m_x0[i];
            dr0_ij[1] = m_y0[j] - m_y0[i];
            dr0_ij[2] = m_z0[j] - m_z0[i];


            for(int d=0; d<m_dim; d++) {
                for(int d2=0; d2<m_dim; d2++) {
                    F(d, d2) += w*dr_ij[d]*dr0_ij[d2]*vol;
                }
            }

            nConnected++;
        }

        if(nConnected <= 3) {
            for(int j=0; j<m_nStressStrainElements; j++) {
                data(i, m_indexStrain[j]) = 0;
                data(i, m_indexStress[j]) = 0;
            }
            continue;
        } else {
            F = F*K; // K = inv(K);
            if(m_greenStrain) {
                strain = 0.5*F.t()*F;
                for(int d=0; d<m_dim;d++) {
                    strain(d, d) -= 0.5;
                }
            } else {
                strain = 0.5*(F.t() + F);
                for(int d=0; d<m_dim;d++) {
                    strain(d, d) -= 1.;
                }
            }

        }

        P(0,0) = (lambda + 2*mu)*strain(0,0) + lambda*strain(1,1) + lambda*strain(2,2);
        P(1,1) = lambda*strain(0,0) + (lambda + 2*mu)*strain(1,1) + lambda*strain(2,2);
        P(2,2) = lambda*strain(0,0) + lambda*strain(1,1) + (lambda + 2*mu)*strain(2,2);
        P(0,1) = 2.*mu*strain(0,1);
        P(0,2) = 2.*mu*strain(0,2);
        P(1,2) = 2.*mu*strain(1,2);

        // Converting PK2 stress to Cauchy stress
        if(m_greenStrain) {
            const double detF = 1./det(F);
            P = detF*F*P*F.t();
        }

        d_stress[0][i] = P(0,0);
        d_stress[1][i] = P(1,1);
        d_stress[2][i] = P(0,1);
        d_stress[3][i] = P(2,2);
        d_stress[4][i] = P(0,2);
        d_stress[5][i] = P(1,2);

        d_strain[0][i] = strain(0,0);
        d_strain[1][i] = strain(1,1);
        d_strain[2][i] = strain(0,1);
        d_strain[3][i] = strain(2,2);
        d_strain[4][i] = strain(0,2);
        d_strain[5][i] = strain(1,2);
#endif
    }
}
//------------------------------------------------------------------------------
void CalculateStressStrain::computeK(int id, int i)
{
//    cout << "Recomputing K " << id << endl;
    const auto &idToCol = m_particles->idToCol();
    const mat &r0 = m_particles->r0();
    mat &data = m_particles->data();
    mat K = zeros(m_dim, m_dim);
    vector<pair<int, vector<double>>> & PDconnections_i = m_particles->pdConnections(id);
    const int nConnections = PDconnections_i.size();
    double dr0_ij[m_dim];

    int nConnected = 0;
    for(int l_j=0; l_j<nConnections; l_j++) {
        const auto &con_i = PDconnections_i[l_j];

//        cout << "id:" << id << " con:" << con_i.first
//             << " connected:" << con_i.second[m_indexConnected]
//                << endl;

        if(con_i.second[m_indexConnected] <= 0.5)
            continue;

        const int id_j = con_i.first;
        const int j = idToCol.at(id_j);

        const double vol_j = data(j, m_indexVolume);
        const double dr0 = con_i.second[m_indexDr0];
        const double volumeScaling_ij = con_i.second[m_indexVolumeScaling];
        const double vol = vol_j*volumeScaling_ij;
        const double w = weightFunction(dr0);

        for(int d=0; d<m_dim; d++) {
            dr0_ij[d] = r0(j, d) - r0(i, d);
        }

        for(int d=0; d<m_dim; d++) {
            for(int d2=0; d2<m_dim; d2++) {
                K(d, d2) += w*dr0_ij[d]*dr0_ij[d2]*vol;
            }
        }

        nConnected++;
    }

//    const vector<int> ids = {56379, 1145, 1208};
//    for(const int id_i:ids) {
//        if(id == id_i) {
//            cout << "id_i: " << id_i << " nConnections:" << nConnections << endl;
//            cout << K << endl;
//        }
//    }

    if(nConnected <= 3) {
        K.zeros();
    } else {
//        K = inv_sympd(K);
//        cout << "before invert id_i: " << id << " nConnections:" << nConnections << "nConnected " << nConnected << endl;
//        cout << K << endl;
        K = inv(K);
    }

//    K(0,1) = 0;
//    K(0,2) = 0;
//    K(1,2) = 0;
//    K(1,0) = 0;
//    K(2,0) = 0;
//    K(2,1) = 0;

    if(m_dim >= 2) {
        data(i, m_indexK[0]) = K(0,0);
        data(i, m_indexK[1]) = K(1,1);
        data(i, m_indexK[2]) = K(0,1);
    }
    if(m_dim == 3) {
        data(i, m_indexK[3]) = K(2,2);
        data(i, m_indexK[4]) = K(0,2);
        data(i, m_indexK[5]) = K(1,2);

//        data(i, m_indexK[0]) = K(0,0);
//        data(i, m_indexK[1]) = K(1,0);
//        data(i, m_indexK[2]) = K(2,0);
//        data(i, m_indexK[3]) = K(0,1);
//        data(i, m_indexK[4]) = K(1,1);
//        data(i, m_indexK[5]) = K(2,1);
//        data(i, m_indexK[6]) = K(0,2);
//        data(i, m_indexK[7]) = K(1,2);
//        data(i, m_indexK[8]) = K(2,2);

    }

//    for(const int id_i:ids) {
//        if(id == id_i) {
//            cout << "id_i: " << id_i << " nConnections:" << nConnections << endl;
//            cout << K << endl;
//        }
//    }
}
//------------------------------------------------------------------------------
}
