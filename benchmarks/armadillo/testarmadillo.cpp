#include <gtest/gtest.h>
#include <PDtools.h>
extern std::vector<string> geometries;
using namespace PDtools;

class PD_ARMADILLO_FIXTURE : public ::testing::Test {
protected:
    PD_Particles testParticles;
    int test_nParticles;
    string loadGeometry;

    PD_ARMADILLO_FIXTURE()
    {
        loadGeometry = geometries[2];
        test_nParticles = 100000;
//        test_nParticles = 60000;
    }
};

//--------------------------------------------------------------------------
// Unodered access
//--------------------------------------------------------------------------
TEST_F(PD_ARMADILLO_FIXTURE, PD_COMPARE_TIMES)
{
    clock_t begin = clock();
    LoadPdParticles loadPdParticles;
    PD_Particles particles = loadPdParticles.load(loadGeometry, "xyz");
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    ASSERT_EQ(particles.nParticles(), test_nParticles);

    double *r_matrix = particles.r().memptr();
    int dim = 3;
    mat &r = particles.r();

    //----------------------------------------------------------------------
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        int i = p_id.second;
        r(0, i) *= 2;
        r(1, i) *=2.1;
        r(2, i) *= 0.1;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Accessing armadillo matrix and changing  all elements took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        double * r_col = r.colptr(p_id.second);
        r_col[0] *= 2;
        r_col[1] *= 2.1;
        r_col[2] *= 0.1;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Using raw pointser to cols took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    // Manually accessing the matrix:
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        int i = p_id.second;
        r_matrix[i*dim] *= 2;
        r_matrix[i*dim + 1] *= 2.1;
        r_matrix[i*dim + 2] *= 0.1;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Accesing and changing the C-matrix elements took " << elapsed_secs << "s" << endl;

    //----------------------------------------------------------------------
    // Manually accessing the matrix:
    begin = clock();
    for(int i=0; i< particles.nParticles(); i++)
    {
        r_matrix[i*dim] *= 2;
        r_matrix[i*dim + 1] *= 2.1;
        r_matrix[i*dim + 2] *= 0.1;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Running over cont. memeory and changing the content of the C-matrix elements took " << elapsed_secs << "s" << endl;

    //----------------------------------------------------------------------
}

//--------------------------------------------------------------------------
// Unodered access
//--------------------------------------------------------------------------
TEST_F(PD_ARMADILLO_FIXTURE, PD_COMPARE_TIMES_RADIUS_UNORDERED)
{
    clock_t begin = clock();
    LoadPdParticles loadPdParticles;
    PD_Particles particles = loadPdParticles.load(loadGeometry, "xyz");
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Loading took " << elapsed_secs << "s" << endl;

    ASSERT_EQ(particles.nParticles(), test_nParticles);

    double *r_matrix = particles.r().memptr();
    int dim = 3;
    mat &r = particles.r();
    arma::ivec A = arma::randi<arma::ivec>(
                particles.nParticles(),
                arma::distr_param(0, particles.nParticles()));

    //----------------------------------------------------------------------
    begin = clock();
    for(int i:A)
    {
        double r_len = sqrt(r(0, i)*r(0, i)
                            + r(1, i)*r(1, i)
                            + r(2, i)*r(2, i));
        r(0, i) = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Axxessing armdillo matrix and changing  all elements took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    begin = clock();
    for(int i:A)
    {
        double * r_col = r.colptr(i);
        double r_len = sqrt(r_col[0]*r_col[0] + r_col[1]*r_col[1] + r_col[2]*r_col[2]);
        r_col[0] = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Using raw pointser to cols took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    // Manually accessing the matrix:
    begin = clock();
    for(int i:A)
    {
        double r_len = sqrt(r_matrix[i*dim]*r_matrix[i*dim]
                + r_matrix[i*dim+1]*r_matrix[i*dim+1]
                + r_matrix[i*dim+2]*r_matrix[i*dim+2]);
        r_matrix[i*dim] = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Accesing and changing the C-matrix elements took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
}

//------------------------------------------------------------------------------
// Ordered access
//------------------------------------------------------------------------------
TEST_F(PD_ARMADILLO_FIXTURE, PD_COMPARE_TIMES_RADIUS_ORDERED)
{
    clock_t begin = clock();
    LoadPdParticles loadPdParticles;
    PD_Particles particles = loadPdParticles.load(loadGeometry, "xyz");
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Loading took " << elapsed_secs << "s" << endl;

    ASSERT_EQ(particles.nParticles(), test_nParticles);

    double *r_matrix = particles.r().memptr();
    int dim = 3;
    mat &r = particles.r();
    //----------------------------------------------------------------------
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        int i = p_id.second;
        double r_len = sqrt(r(0, i)*r(0, i)
                            + r(1, i)*r(1, i)
                            + r(2, i)*r(2, i));
        r(0, i) = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Axxessing armdillo matrix and changing  all elements took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        int i = p_id.second;
        double * r_col = r.colptr(i);
        double r_len = sqrt(r_col[0]*r_col[0] + r_col[1]*r_col[1] + r_col[2]*r_col[2]);
        r_col[0] = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Using raw pointser to cols took " << elapsed_secs << "s" << endl;
    //----------------------------------------------------------------------
    // Manually accessing the matrix:
    begin = clock();
    for(auto p_id:particles.pIds())
    {
        int i = p_id.second;
        double r_len = sqrt(r_matrix[i*dim]*r_matrix[i*dim]
                + r_matrix[i*dim+1]*r_matrix[i*dim+1]
                + r_matrix[i*dim+2]*r_matrix[i*dim+2]);
        r_matrix[i*dim] = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Accesing and changing the C-matrix elements took " << elapsed_secs << "s" << endl;

    //----------------------------------------------------------------------
    // Manually accessing the matrix:
    begin = clock();
    for(int i=0; i< particles.nParticles(); i++)
    {
        double r_len = sqrt(r_matrix[i*dim]*r_matrix[i*dim]
                + r_matrix[i*dim+1]*r_matrix[i*dim+1]
                + r_matrix[i*dim+2]*r_matrix[i*dim+2]);
        r_matrix[i*dim] = r_len;
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Running over cont. memeory and changing the content of the C-matrix elements took " << elapsed_secs << "s" << endl;
}
