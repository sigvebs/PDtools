
#include <gtest/gtest.h>
#include <PdFunctions/pdfunctions.h>
#include <PDtools.h>
#include <PDtools/Solver/solvers.h>
#include <PDtools/Force/forces.h>

#include <PDtools/Modfiers/modifiers.h>

extern std::vector<string> geometries;

using namespace PDtools;

class PD_LINEAR_SOLVER_FIXTURE : public ::testing::Test {
protected:

    PD_LINEAR_SOLVER_FIXTURE()
    {

    }
};

TEST_F(PD_LINEAR_SOLVER_FIXTURE, TEST_CG)
{
    using namespace std;
    using namespace arma;

    int m_dim = 2;

    const double appliedForce = 2e-17;
    const double threshold = 3.e-8;
    const double vMag = 1e-11;
    const double s0 = 0.00158377370467;
    const double rho = 2140.0;
    double lc = 2.66216165e-06/3;
    const double delta = 2.66216165e-06;
    const double gridspacing = 1.1*delta;

    const double h = 1.25e-05;
    const double E =  1.52e10;
    const double nu = 1./3.;
    const double k = E/(2.*(1. - nu));
    double micromodulus = 12.*k/(M_PI*h*pow(delta, 3));
    micromodulus = k;
    vector<string> saveParameters = {"id", "x", "y", "potential_energy", "kinetic_energy", "stress", "damage"};

    const int maxIterations = 20000;
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------
    PD_Particles m_particles;
    m_particles = load_pd(geometries[8]);
    m_particles.initializeADR();
    m_particles.registerParameter("rho", rho);
    const int m_indexMicromodulus = m_particles.registerParameter("micromodulus", micromodulus);
    m_particles.registerParameter("s0", s0);
    m_particles.registerParameter("radius", s0);

    const int m_nParticles = m_particles.nParticles();

    // New variables needed for solving the lienar equations us CG
    const int degFreedom = m_dim*m_nParticles;
    arma::vec u_k = arma::vec(degFreedom);
    arma::vec r_k = arma::vec(degFreedom);
    arma::vec r_k1 = arma::vec(degFreedom);
    arma::vec p_k = arma::vec(degFreedom);

    vector<pair<double,double>> boundaries;

    pair<double,double> x_limits(0., 5e-05);
    pair<double,double> y_limits(0., 5e-05);
    pair<double,double> z_limits(-0.5*h, 0.5*h);

    boundaries = {x_limits, y_limits, z_limits};

    Grid grid(boundaries, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(m_particles);
    calculateRadius(m_particles, m_dim, h);

    lc *= 1.05;
    setPdConnections(m_particles, grid, delta, lc);
    m_particles.registerPdParameter("volumeScaling", 1);
    applyVolumeCorrection(m_particles, delta, lc);
    //--------------------------------------------------------------------------
    // Setting the initial position
    //--------------------------------------------------------------------------
    arma::mat &m_dr0 = m_particles.r0();
    arma::mat &m_r = m_particles.r();

    for(int i=0; i<m_particles.nParticles(); i++)
    {
        for(int d=0; d<3; d++)
        {
            m_dr0(d, i) = m_r(d,i);
        }
    }
    //--------------------------------------------------------------------------
    // Setting the stiffness matrix
    //--------------------------------------------------------------------------
    arma::mat & m_data = m_particles.data();
    std::unordered_map<int, int> & m_pIds = m_particles.pIds();

    const int m_indexVolume = m_particles.getParamId("volume");
    const int m_indexDr0 = m_particles.getPdParamId("dr0");
    const int m_indexVolumeScaling = m_particles.getPdParamId("volumeScaling");
    const int m_indexConnected = m_particles.getPdParamId("connected");

    cout << "Building stiffness matrix" << endl;
#if 0
    arma::sp_mat C(degFreedom, degFreedom);
    for(int l_a=0; l_a<m_nParticles; l_a++)
    {
        pair<int, int> id(l_a, l_a);
        const int pId = id.first;
        const int a = id.second;
        const double c_i = m_data(a, m_indexMicromodulus);

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        const int nConnections = PDconnections.size();

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            const auto &con = PDconnections[l_j];
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int b = m_pIds[id_j];

            const double c_j = m_data(b, m_indexMicromodulus);
            const double vol_j = m_data(b, m_indexVolume);
            const double dr0         = con.second[m_indexDr0];
            const double volumeScaling   = con.second[m_indexVolumeScaling];
            const double c_ab = 0.5*(c_i + c_j);

            const double coeff = c_ab/(pow(dr0, 3))*vol_j*volumeScaling;
            arma::vec dr0_v = m_dr0.col(l_a) - m_dr0.col(b);

            for(int d1=0; d1<m_dim; d1++)
            {
                int i_a = a*m_dim;
                int i_b = b*m_dim;

                // Alpha-alpha sum
                C(i_a + d1, i_a + d1) += coeff*dr0_v(d1)*dr0_v(d1);

                // Alpha-beta
                C(i_a + d1, i_b + d1) -= coeff*dr0_v(d1)*dr0_v(d1);

                for(int d2=d1+1; d2<m_dim; d2++)
                {
                    double element = coeff*dr0_v(d1)*dr0_v(d2);

                    // Alpha-alpha sum
                    C(i_a + d1, i_a + d2) += element;
                    C(i_a + d2, i_a + d1) += element;

                    // Alpha-beta
                    C(i_a + d1, i_b + d2) -= element;
                    C(i_a + d2, i_b + d1) -= element;
                }
            }
        }
    }
#else
    // Finding the size of locations and values

    int total_values = 0;
    for(int l_a=0; l_a<m_nParticles; l_a++)
    {
        pair<int, int> id(l_a, l_a);
        const int pId = id.first;
        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        const int nConnections = PDconnections.size();
        total_values += nConnections + 1;
    }

    arma::umat locations(2, total_values*m_dim*m_dim);
    arma::vec values(total_values*m_dim*m_dim);
    int pos = 0;


    for(int l_a=0; l_a<m_nParticles; l_a++)
    {
        pair<int, int> id(l_a, l_a);
        const int pId = id.first;
        const int a = id.second;
        int i_a = a*m_dim;
        const double c_i = m_data(a, m_indexMicromodulus);

        vector<pair<int, vector<double>>> & PDconnections = m_particles.pdConnections(pId);
        const int nConnections = PDconnections.size();
        arma::mat C_ij = arma::zeros(m_dim, m_dim);

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            const auto &con = PDconnections[l_j];
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int b = m_pIds[id_j];
            const int i_b = b*m_dim;

            const double c_j = m_data(b, m_indexMicromodulus);
            const double vol_j = m_data(b, m_indexVolume);
            const double dr0         = con.second[m_indexDr0];
            const double volumeScaling   = con.second[m_indexVolumeScaling];
            const double c_ab = 0.5*(c_i + c_j);

            const double coeff = c_ab/(pow(dr0, 3))*vol_j*volumeScaling;
            arma::vec dr0_v = m_dr0.col(l_a) - m_dr0.col(b);

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

    arma::sp_mat C(locations, values, degFreedom, degFreedom);
#endif
    cout << "Stiffness matrix complete" << endl;
    //--------------------------------------------------------------------------
    // Setting the boundary conditions
    //--------------------------------------------------------------------------
    arma::vec b = arma::zeros(degFreedom);
    pair<double, double> m_boundary_left(x_limits.first - delta, x_limits.first + delta);
    pair<double, double> m_boundary_right(x_limits.second - delta, x_limits.second + delta);
    int m_boundaryOrientation = 0;
    double appliedVolume = h*(m_boundary_left.second - m_boundary_left.first);
    double bodyForce = appliedForce/appliedVolume;

    for(int i=0; i<m_particles.nParticles(); i++)
    {
        int col_i = i;
        double pos = m_r(m_boundaryOrientation, col_i);
        if(m_boundary_left.first <= pos && pos < m_boundary_left.second)
        {
            pair<int, int> pId(i, i);
            int j = i*m_dim + m_boundaryOrientation;
            b(j) = -bodyForce;
        }
        if(m_boundary_right.first <= pos && pos < m_boundary_right.second)
        {
            pair<int, int> pId(i, i);
            int j = i*m_dim + m_boundaryOrientation;
            b(j) = bodyForce;
        }
    }
    //--------------------------------------------------------------------------
    // Setting the initial displacement
    //--------------------------------------------------------------------------
    u_k.zeros();
    double stretch = 0.01;

    for(int i=0; i<m_particles.nParticles(); i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            int j = i*m_dim + d;
            u_k(j) = stretch*m_dr0(d, i);
        }
    }
//    u_k.randu();
//    u_k *= 0.75*lc;
    //--------------------------------------------------------------------------
    // CG algorithm
    //--------------------------------------------------------------------------
    r_k = b - C*u_k;
    p_k = r_k;
    arma::vec Ap(degFreedom);
    int l_k;
    double r_k2;

    for(l_k=0; l_k<maxIterations; l_k++)
    {
        Ap = C*p_k;
        r_k2 = arma::dot(r_k, r_k);
        const double alpha = r_k2/(arma::dot(p_k, Ap));
        r_k1 = r_k - alpha*Ap;
        u_k += alpha*p_k;
        const double beta = arma::dot(r_k1, r_k)/r_k2;
        if(sqrt(r_k2) < threshold)
             break;
        p_k = r_k1 + beta*p_k;
        r_k = r_k1;
//        arma::swap(r_k1, r_k);
    }
    cout << "Error: " << sqrt(r_k2) << " at iteration " << l_k << endl;

    //--------------------------------------------------------------------------
    // Setting the final position
    //--------------------------------------------------------------------------
//    arma::mat u_vec(u_k.memptr(), m_dim, m_nParticles);
    for(int i=0; i<m_particles.nParticles(); i++)
    {
        for(int d=0; d<m_dim; d++)
        {
            int j = i*m_dim + d;
            m_r(d,i) = m_dr0(d, i) + u_k(j);
        }
    }

    //--------------------------------------------------------------------------
    // Saving the results
    //--------------------------------------------------------------------------
    string save_folder = "/media/Media4/Scratch/Untitled Folder";
    save_folder = "/media/sigve/Pengebingen/scratch/tmp";
    SaveParticles save_lmp("lmp");
    string save_path = save_folder + "/test" + to_string(m_nParticles) + ".lmp";
    save_lmp.writeToFile(m_particles, save_path);
    /*
    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    delete solver;
    delete pdForce;
    delete boundaryLeft;
    delete boundaryRight;
    */
}
