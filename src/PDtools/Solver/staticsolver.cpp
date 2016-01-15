#include "staticsolver.h"

#include "PDtools/Particles/pd_particles.h"
#include "PDtools/Force/force.h"
#include "PDtools/Modfiers/modifier.h"
#include "PDtools/SavePdData/savepddata.h"


namespace PDtools
//------------------------------------------------------------------------------
{
StaticSolver::StaticSolver(int maxIterations, int threshold):
    m_maxIterations(maxIterations),
    m_threshold(threshold)
{

}
//------------------------------------------------------------------------------
void StaticSolver::solve()
{
    initialize();
    checkInitialization();

    // Looping over all time, particles and components.
    for (int i = 0; i < m_steps; i++)
    {
        stepForward(i);
    }
}
//------------------------------------------------------------------------------
void StaticSolver::stepForward(int i)
{
    save(i);
    modifiersStepOne();

    iterate();
//    updateStiffnessMatrix();

    cout << "i = " << i << endl;

    modifiersStepTwo();
    m_t += m_dt;
}
//------------------------------------------------------------------------------
void StaticSolver::iterate()
{
//ARMA_USE_WRAPPER
#ifdef ARMA_USE_WRAPPER
    const int nParticles = m_particles->nParticles();
    //--------------------------------------------------------------------------
    // Setting the boundary conditions
    //--------------------------------------------------------------------------
//    vec b = conv_to<vec>::from(m_particles->b()); // TMP!
    const mat &B = m_particles->b();
    mat &U = m_particles->u();
//    B.reshape(nParticles*m_dim, 1);
//    vec b_k = B.col(0);

    for(int i=0; i<nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            const int j = i*m_dim + d;
            b_k(j) = B(i, d);
            u_k(j) = U(i, d);
        }
    }
    //--------------------------------------------------------------------------
    // CG
    //--------------------------------------------------------------------------
    vec test = C*u_k;
    r_k = b_k - C*u_k;
    p_k = r_k;
    int l_k = 0;
    double r_k2;
    r_k2 = dot(r_k, r_k);
    double t2 = r_k2 = dot(test, test);
    cout << "Error1: " << sqrt(r_k2) << " at iteration " << l_k << endl;
    cout << "Error2: " << sqrt(t2) << " at iteration " << l_k << endl;
    cout << "Error3: " << sqrt(dot(b_k, b_k)) << " at iteration " << l_k << endl;

    m_maxIterations = 2000;
    for(l_k=0; l_k<m_maxIterations; l_k++)
    {
        Ap = C*p_k;
        r_k2 = dot(r_k, r_k);
        const double alpha = r_k2/(arma::dot(p_k, Ap));
        r_k1 = r_k - alpha*Ap;
        u_k += alpha*p_k;
        const double beta = arma::dot(r_k1, r_k)/r_k2;
        if(sqrt(r_k2) < m_threshold)
             break;
        p_k = r_k1 + beta*p_k;
//        r_k = r_k1;
        arma::swap(r_k1, r_k);
    }
    r_k2 = dot(r_k, r_k);
    cout << "Error: " << sqrt(r_k2) << " at iteration " << l_k << endl;

    //--------------------------------------------------------------------------
    // Setting the final position
    //--------------------------------------------------------------------------
    mat &r = m_particles->r();
    mat &r0 = m_particles->r0();

    for(unsigned int i=0; i<nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            const int j = i*m_dim + d;
            U(i, d) = u_k(j);
            r(i, d) = r0(i, d) + U(i, d);
        }
    }
#endif
}
//------------------------------------------------------------------------------
void StaticSolver::initialize()
{
    const int nParticles = m_particles->nParticles();
    m_particles->initializeBodyForces();
    m_degFreedom = m_dim*nParticles;
    u_k = arma::vec(m_degFreedom);
    r_k = arma::vec(m_degFreedom);
    r_k1 = arma::vec(m_degFreedom);
    p_k = arma::vec(m_degFreedom);
    b_k = arma::vec(m_degFreedom);
    Ap = arma::vec(m_degFreedom);

    //--------------------------------------------------------------------------
    // Setting the initial displacement
    //--------------------------------------------------------------------------
    mat & u = m_particles->u();
    const double stretch = 1e-9;
    const mat &r0 = m_particles->r0();
    for(int i=0; i<nParticles; i++)
    {
        for(int d=0; d<m_dim; d++)
        {
//            const int j = i*m_dim + d;
//            u_k(j) = stretch*r0(i, d);
            u(i, d) = stretch*r0(i, d);
            u(i, d) = 0;
        }
    }

    //-----------------------------
    createStiffnessMatrix();
    arma::mat & b = m_particles->b();
    b.zeros();
    Solver::initialize();
}
//------------------------------------------------------------------------------
void StaticSolver::save(int i)
{
    if(i%m_saveInterval == 0)
    {
        m_saveParticles->evaluate(m_t, i);
        computeStress();
        m_saveParticles->saveData(m_t, i);
    }
}
//------------------------------------------------------------------------------
void StaticSolver::createStiffnessMatrix()
{
    const int nParticles = m_particles->nParticles();
    arma::mat & m_data = m_particles->data();
    std::unordered_map<int, int> & m_idToCol = m_particles->idToCol();

    arma::mat &m_dr0 = m_particles->r0();

    const int m_indexVolume = m_particles->getParamId("volume");
    const int m_indexDr0 = m_particles->getPdParamId("dr0");
    const int m_indexVolumeScaling = m_particles->getPdParamId("volumeScaling");
    const int m_indexConnected = m_particles->getPdParamId("connected");
    const int m_indexMicromodulus = m_particles->getParamId("micromodulus");

    int total_values = 0;
    for(int i=0; i<nParticles; i++)
    {
        pair<int, int> id(i, i);
        const int pId = id.first;
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
        const int nConnections = PDconnections.size();
        total_values += nConnections + 1;
    }

    arma::umat locations(2, total_values*m_dim*m_dim);
    arma::vec values(total_values*m_dim*m_dim);
    int pos = 0;


    for(int l_a=0; l_a<nParticles; l_a++)
    {
        pair<int, int> id(l_a, l_a);
        const int pId = id.first;
        const int a = id.second;
        int i_a = a*m_dim;
        const double c_i = m_data(a, m_indexMicromodulus);

        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(pId);
        const int nConnections = PDconnections.size();
        arma::mat C_ij = arma::zeros(m_dim, m_dim);

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            const auto &con = PDconnections[l_j];
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int b = m_idToCol.at(id_j);
            const int i_b = b*m_dim;

            const double c_j = m_data(b, m_indexMicromodulus);
            const double vol_j = m_data(b, m_indexVolume);
            const double dr0         = con.second[m_indexDr0];
            const double volumeScaling   = con.second[m_indexVolumeScaling];
            const double c_ab = 0.5*(c_i + c_j);

            const double coeff = c_ab/(pow(dr0, 3))*vol_j*volumeScaling;
            const vec & r_la = m_dr0.row(l_a).t();
            const vec & r_b = m_dr0.row(b).t();
            const arma::vec dr0_v = r_la - r_b;

            for(int d1=0; d1<m_dim; d1++)
            {
                // Alpha-alpha sum
                C_ij(d1, d1) += coeff*dr0_v(d1)*dr0_v(d1);

                // Alpha-beta
                locations(0, pos) = i_a + d1;
                locations(1, pos) = i_b + d1;
                values(pos++) = -coeff*dr0_v(d1)*dr0_v(d1);

                for(int d2=d1+1; d2<m_dim; d2++)
                {
                    const double element = coeff*dr0_v(d1)*dr0_v(d2);

                    // Alpha-alpha sum
                    C_ij(d1, d2) += element;
                    C_ij(d2, d1) += element;

                    // Alpha-beta
                    locations(0, pos) = i_a + d1;
                    locations(1, pos) = i_b + d2;
                    values(pos++) = -element;
                    locations(0, pos) = i_a + d2;
                    locations(1, pos) = i_b + d1;
                    values(pos++) = -element;
                }
            }
        }
        for(int d1=0; d1<m_dim; d1++)
        {
            locations(0, pos) = i_a + d1;
            locations(1, pos) = i_a + d1;
            values(pos++) = C_ij(d1, d1);

            for(int d2=d1+1; d2<m_dim; d2++)
            {
                locations(0, pos) = i_a + d1;
                locations(1, pos) = i_a + d2;
                values(pos++) = C_ij(d1, d2);

                locations(0, pos) = i_a + d2;
                locations(1, pos) = i_a + d1;
                values(pos++) = C_ij(d2, d1);
            }
        }
    }

    C = arma::sp_mat(locations, values, m_degFreedom, m_degFreedom);
    cout << "Stiffness matrix complete" << endl;
}
//------------------------------------------------------------------------------
void StaticSolver::updateStiffnessMatrix()
{

}
//------------------------------------------------------------------------------
void StaticSolver::computeStress()
{
    int indexStress[6];
    const int m_indexCompute = m_particles->getPdParamId("compute");
    const int m_indexVolume = m_particles->getParamId("volume");
    const int m_indexDr0 = m_particles->getPdParamId("dr0");
    const int m_indexVolumeScaling = m_particles->getPdParamId("volumeScaling");
    const int m_indexConnected = m_particles->getPdParamId("connected");
    const int m_indexMicromodulus = m_particles->getParamId("micromodulus");
    const int m_indexForceScaling = m_particles->getPdParamId("forceScalingBond");

    indexStress[0] = m_particles->getParamId("s_xx");
    indexStress[1] = m_particles->getParamId("s_yy");
    indexStress[2] = m_particles->getParamId("s_zz");
    indexStress[3] = m_particles->getParamId("s_xy");
    indexStress[4] = m_particles->getParamId("s_xz");
    indexStress[5] = m_particles->getParamId("s_yz");
    arma::mat & m_data = m_particles->data();
    const arma::mat & U = m_particles->u();
    const arma::mat & R = m_particles->r();
    const arma::mat &m_dr0 = m_particles->r0();
    const std::unordered_map<int, int> & idToCol= m_particles->idToCol();

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(unsigned int a=0; a<m_particles->nParticles(); a++)
    {
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(a);
        const int nConnections = PDconnections.size();
        const double c_i = m_data(a, m_indexMicromodulus);
        arma::mat C_ij = arma::zeros(m_dim, m_dim);

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            const auto &con = PDconnections.at(l_j);
            if(con.second[m_indexConnected] <= 0.5 || !con.second[m_indexCompute])
                continue;
            const int id_b = con.first;
            const int b = idToCol.at(id_b);

            const double c_j = m_data(b, m_indexMicromodulus);
            const double vol_j = m_data(b, m_indexVolume);
            const double dr0         = con.second[m_indexDr0];
            const double volumeScaling   = con.second[m_indexVolumeScaling];
            const double g_ij = con.second[m_indexForceScaling];
            const double c_ab = 0.5*(c_i + c_j)*g_ij;

            const double coeff = c_ab/(pow(dr0, 3))*vol_j*volumeScaling;
            const vec & r_la = m_dr0.row(a).t();
            const vec & r_b = m_dr0.row(b).t();
            const arma::vec dr0_v = r_la - r_b;
            for(int d1=0; d1<m_dim; d1++)
            {
                // Alpha-alpha sum
                C_ij(d1, d1) = coeff*dr0_v(d1)*dr0_v(d1);

                for(int d2=d1+1; d2<m_dim; d2++)
                {
                    const double element = coeff*dr0_v(d1)*dr0_v(d2);
                    C_ij(d1, d2) = element;
                    C_ij(d2, d1) = element;
                }
            }

            vec u_ab = U.row(b).t() - U.row(a).t();
            vec dr_ab = R.row(b).t() - R.row(a).t();
            vec k(m_dim);
            k.zeros();

            for(int d1=0;d1<m_dim; d1++)
            {
                for(int d2=0;d2<m_dim; d2++)
                {
                    k(d1) += C_ij(d1, d2)*u_ab(d2);
                }
            }

            m_data(a, indexStress[0]) += 0.5*k[0]*dr_ab[0];
            m_data(a, indexStress[1]) += 0.5*k[1]*dr_ab[1];
            m_data(a, indexStress[3]) += 0.5*k[0]*dr_ab[1];

            if(m_dim == 3)
            {
                m_data(a, indexStress[2]) += 0.5*k[2]*dr_ab[2];
                m_data(a, indexStress[4]) += 0.5*k[0]*dr_ab[2];
                m_data(a, indexStress[5]) += 0.5*k[1]*dr_ab[2];
            }
        }
    }
}
//------------------------------------------------------------------------------
}
