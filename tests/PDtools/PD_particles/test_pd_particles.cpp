#include <gtest/gtest.h>
#include <PDtools.h>
#include <PdFunctions/pdfunctions.h>

extern std::vector<string> geometries;

using namespace PDtools;

class PD_PARTICLES_FIXTURE : public ::testing::Test {
    protected:
    PD_Particles testParticles;
    int test_nParticles = 11;

    PD_PARTICLES_FIXTURE()
    {
        testParticles.nParticles(test_nParticles);
        testParticles.initializeMatrices();
        mat & data = testParticles.data();
        mat & r = testParticles.r();
        mat & v = testParticles.v();

        auto & parameters = testParticles.parameters();
        parameters["volume"] = 0;
        arma::vec volume;

        mat _r, _v;

        _r << 0.313150 << 0.932484 << 0.731905 << 0.650592 << 0.109496 << 0.532164 << 0.685139 << 0.678746 << 0.384222 << 0.947418 << 0.294880 << arma::endr
                        << 0.597834 << 0.519061 << 0.547958 << 0.343276 << 0.449608 << 0.704623 << 0.468903 << 0.390759 << 0.581632 << 0.401296 << 0.368566 << arma::endr;
        _v << 0.165347 << 0.001397  << 0.212907 << 0.529744  << 0.195344 << 0.504753 << 0.236449 << 0.546610 << 0.681405 << 0.665064 << 0.110014 << arma::endr
                        << 0.265347 << 0.101397 << 0.112907 << 0.329744 << 0.695344 << 0.0504753 << 0.536449 << 0.346610 << 0.168145 << 0.166564 << 0.111014 << arma::endr;
        volume << 8.76566e-06 << 8.22899e-06 << 8.0501e-06  << 8.31843e-06 << 7.78176e-06 << 6.44008e-06 << 7.33453e-06 << 8.0501e-06 << 8.31843e-06 << 7.0662e-06 << 8.13954e-06 << arma::endr;

        r = _r.t();
        v = _v.t();
        data.col(parameters["volume"]) = volume;

        unordered_map<int, int> & pIds = testParticles.pIds();
        for(int i=0; i<test_nParticles;i++)
        {
            pIds[i] = i;
        }
    }
};


TEST_F(PD_PARTICLES_FIXTURE, PD_LOAD_XYZ)
{
    LoadPdParticles loadPdParticles;
    PD_Particles particles = loadPdParticles.load(geometries[0], "xyz");

    ASSERT_EQ(particles.nParticles(), 11)
            << "The number of particles is wrong.";;
    ASSERT_EQ(particles.parameters().count("volume"), 1)
            << "The 'volume' parameter is not set for the particles.";;

    const mat & r = particles.r();
    const mat & r_test = testParticles.r();
    const mat & v = particles.v();
    const mat & v_test = testParticles.v();
    const mat & data = particles.data();
    const mat & data_test = testParticles.data();
    const auto & parameters = particles.parameters();
    auto & parameters_test = testParticles.parameters();

    for(int i=0; i<test_nParticles; i++)
    {
        int id_test =  testParticles.pIds()[i];
        int id_load =  particles.pIds()[i];

        ASSERT_EQ(r_test(id_test, 0), r(id_load, 0));
        ASSERT_EQ(r_test(id_test, 1), r(id_load, 1));
        ASSERT_EQ(v_test(id_test, 0), v(id_load, 0));
        ASSERT_EQ(v_test(id_test, 1), v(id_load, 1));

        for(auto param:parameters)
        {
            int particles_type_pos = param.second;
            int testParticles_type_pos = parameters_test[param.first];
            ASSERT_EQ(data_test(id_test, testParticles_type_pos), data(id_load, particles_type_pos))
                    << "The data in Testparticles does not match the data read from file";
        }
    }

}

TEST_F(PD_PARTICLES_FIXTURE, PD_SAVE_AND_LOAD_XYZ)
{
    LoadPdParticles loadPDParticles;
    PD_Particles particles = loadPDParticles.load(geometries[2], "xyz");



    SaveParticles save_lmp("lmp");
    save_lmp.writeToFile(particles, string(TEST_SAVE_PATH) + "/testFromPD.lmp");
}

TEST_F(PD_PARTICLES_FIXTURE, CREATE_PD_CONNECTIONS)
{
    using namespace std;
    using namespace arma;

    LoadPdParticles loadParticles;
    PD_Particles particles;
    particles = loadParticles.load(geometries[0], "xyz");
//    particles = loadParticles.load(geometries[4], "xyz");

    vector<pair<double,double>> boundaries;

    // Dots - geometry[0]
    pair<double,double> x_limits(0., 1.);
    pair<double,double> y_limits(0., 0.9);
    pair<double,double> z_limits(0., 1.);
    double gridspacing = 0.18*2;
    double radius = 0.95*gridspacing;

    // Bunny - geometry[1]
    //    pair<double,double> x_limits(0., 1.);
    //    pair<double,double> y_limits(0., 0.9);
    //    pair<double,double> z_limits(0., 1.);
    //    double gridspacing = 0.018*2;
    //    double radius = 0.95*gridspacing;

    //    // funny - geometry[2]
    //    pair<double,double> x_limits(0., 1.);
    //    pair<double,double> y_limits(0., 0.9);
    //    pair<double,double> z_limits(0., 1.);
    //    double gridspacing = 0.008*2.5;
    //    double radius = 0.95*gridspacing;

    // Hole - geometry[4]
    //    pair<double,double> x_limits(0., 0.0487619);
    //    pair<double,double> y_limits(0., 0.0488042);
    //    pair<double,double> z_limits(0., 0.001);
    //    double gridspacing = 0.0022*3;;
    //    double radius = 0.95*gridspacing;

    boundaries = {x_limits, y_limits, z_limits};

    Grid grid(boundaries, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(particles);

    setPdConnections(particles, grid, radius, gridspacing);

    // Saving the results
    particles.registerParameter("connections");
    int connectionsID = particles.parameters("connections");
    particles.registerParameter("gridID");
    int gridID = particles.parameters("gridID");

    const unordered_map<int, int> & pIds = particles.pIds();
    const mat & R = particles.r();
    mat & data = particles.data();

    for(const pair<int, int> &idCol:pIds)
    {
        int gridId = grid.gridId(R.row(idCol.second));
        data(idCol.second, gridID) = gridId;
        data(idCol.second, connectionsID) = particles.pdConnections(idCol.first).size();
    }
    SaveParticles save_lmp("lmp");
    save_lmp.writeToFile(particles, string(TEST_SAVE_PATH) + "/pdConnections.lmp");
    SaveParticles save_xyz("xyz");
    save_xyz.writeToFile(particles, string(TEST_SAVE_PATH) + "/pdConnections.xyz");
}
