#include <gtest/gtest.h>
#include <PDtools.h>

using namespace PDtools;
using namespace arma;

extern std::vector<string> geometries;

class PARTICLES_FIXTURE : public ::testing::Test {
protected:
    Particles testParticles;
    int test_nParticles = 11;

    PARTICLES_FIXTURE()
    {
        testParticles.nParticles(test_nParticles);
        testParticles.initializeMatrices();
        mat & data = testParticles.data();
        mat & r = testParticles.r();

        auto & parameters = testParticles.parameters();
        parameters["v_x"] = 0;
        parameters["v_y"] = 1;
        parameters["volume"] = 2;

        vec v_x, v_y, volume;
        mat _r;
        _r << 0.313150 << 0.932484 << 0.731905 << 0.650592 << 0.109496 << 0.532164 << 0.685139 << 0.678746 << 0.384222 << 0.947418 << 0.294880 << arma::endr
                        << 0.597834 << 0.519061 << 0.547958 << 0.343276 << 0.449608 << 0.704623 << 0.468903 << 0.390759 << 0.581632 << 0.401296 << 0.368566 << arma::endr;
        v_x << 0.165347 << 0.001397  << 0.212907 << 0.529744  << 0.195344 << 0.504753 << 0.236449 << 0.546610 << 0.681405 << 0.665064 << 0.110014;
        v_y << 0.265347 << 0.101397 << 0.112907 << 0.329744 << 0.695344 << 0.0504753 << 0.536449 << 0.346610 << 0.168145 << 0.166564 << 0.111014;
        volume << 8.76566e-06 << 8.22899e-06 << 8.0501e-06  << 8.31843e-06 << 7.78176e-06 << 6.44008e-06 << 7.33453e-06 << 8.0501e-06 << 8.31843e-06 << 7.0662e-06 << 8.13954e-06 << arma::endr;

        r = _r;
        data.col(parameters["v_x"]) = v_x;
        data.col(parameters["v_y"]) = v_y;
        data.col(parameters["volume"]) = volume;

        unordered_map<int, int> & pIds = testParticles.pIds();
        for(int i=0; i<test_nParticles;i++)
        {
            pIds[i] = i;
        }
    }
};

TEST_F(PARTICLES_FIXTURE, LOAD_XYZ)
{
    LoadParticles loadParticles;
    Particles particles = loadParticles.load(geometries[0], "xyz");
    cout << "loaded" << endl;

    ASSERT_EQ(particles.nParticles(), test_nParticles)
            << "The number of particles is wrong.";

    ASSERT_EQ(particles.parameters().count("volume"), 1)
            << "The 'volume' parameter is not set for the particles.";;

    const mat & r = particles.r();
    const mat & r_test = testParticles.r();
    const mat & data = particles.data();
    const mat & data_test = testParticles.data();
    const auto & parameters = particles.parameters();
    auto & parameters_test = testParticles.parameters();


    for(int i=0; i<test_nParticles; i++)
    {
        int id_test = testParticles.pIds()[i];
        int id_load = particles.pIds()[i];

        ASSERT_EQ(r_test(0, id_test), r(0, id_load));
        ASSERT_EQ(r_test(1, id_test), r(1, id_load));
        ASSERT_EQ(r_test(0, id_test), r(0, id_load));
        ASSERT_EQ(r_test(1, id_test), r(1, id_load));

        for(auto param:parameters)
        {
            int particles_type_pos = param.second;
            int testParticles_type_pos = parameters_test[param.first];

            ASSERT_EQ(data_test(id_test, testParticles_type_pos), data(id_load, particles_type_pos))
                    << "The data in Testparticles does not match the data read from file";
        }
    }
}

TEST_F(PARTICLES_FIXTURE, LOAD_PLY)
{
    LoadParticles loadParticles;
//    Particles particles = loadParticles.load(geometries[7], "ply", true);
}

TEST_F(PARTICLES_FIXTURE, PLACE_PARTICLES_IN_GRID)
{
    using namespace std;
    using namespace arma;

    vector<pair<double,double>> boundaries;
    pair<double,double> x_limits(0., 1.);
    pair<double,double> y_limits(0., .8);
    pair<double,double> z_limits(0., 1.);
    boundaries = {x_limits, y_limits, z_limits};
    double gridspacing = 0.016*3;
//    double gridspacing = 0.1;

    Grid grid(boundaries, gridspacing);
    grid.initialize();


    //    LoadParticles loadParticles;
    unordered_map<string, int> inputParameter =
                                    {{"x",0}, {"y",1,}, {"z",2},
                                     {"s_xz", 3}, {"s_xy", 4},
                                     {"s_y", 5}, {"s_x", 6},
                                     {"s_z", 7},
                                     {"damage", 8},
                                     {"grid_id", 9},
                                     {"stress_I", 10},
                                     {"kinetic_energy", 11},
                                     {"F_x", 12},
                                     {"id", 13},
                                     {"s_yz", 14}
                                    };


    //    Particles particles = loadParticles.load(geometries[2], "xyz");
    LoadParticles loadParticles;
    loadParticles.useLegacyFormat(true);
    Particles particles;
//    particles = loadParticles.load(geometries[6], "lmp", true, inputParameter);
    particles = loadParticles.load(geometries[1], "xyz");


//    grid.placeParticlesInGrid(particles);
////    vector<string> saveParameters = {"x", "y", "z"};
//    vector<string> saveParameters = {"x", "y", "z",  "volume"};

//    SaveParticles save_xyz("xyz", saveParameters);
//    save_xyz.writeToFile(particles, string(TEST_SAVE_PATH) + "/test.xyz");
//    SaveParticles save_lmp("lmp");
//    save_lmp.writeToFile(particles, string(TEST_SAVE_PATH) + "/test.lmp");
//    SaveParticles save_ply("ply");
//    save_ply.writeToFile(particles, string(TEST_SAVE_PATH) + "/test.ply");

//    SaveParticles save_lmpBinary("lmp", true);
//    save_lmpBinary.writeToFile(particles, string(TEST_SAVE_PATH) + "/testBinary.dump");
//    SaveParticles save_plyBinary("ply", true);
////    save_plyBinary.append(true);
//    save_plyBinary.writeToFile(particles, string(TEST_SAVE_PATH) + "/testBinary.ply");

//    particles = loadParticles.load(string(TEST_SAVE_PATH) + "/testBinary.ply", "ply", true);
//    SaveParticles test("lmp");
//    test.writeToFile(particles, string(TEST_SAVE_PATH) + "/testFromPly.lmp");

//    //    unordered_map<string, int> inputParameters = {"id", "x", "y", "z", "volume"};
//    inputParameter = {{"x", 1}, {"y", 2}, {"z", 3}, {"id", 0}, {"damage", 4}};
//    loadParticles.useLegacyFormat(false);
//    particles = loadParticles.load(string(TEST_SAVE_PATH) + "/testBinary.dump", "lmp", true, inputParameter);
//    save_lmp.writeToFile(particles, string(TEST_SAVE_PATH) + "/testFromBin.lmp");
}


TEST_F(PARTICLES_FIXTURE, VERLET_LIST)
{
    using namespace std;
    using namespace arma;

    LoadParticles loadParticles;
    Particles particles;
    particles = loadParticles.load(geometries[2], "xyz");
//    particles = loadParticles.load(geometries[4], "xyz");
    string verletStringId = "verletList";
    particles.registerVerletList(verletStringId);

    vector<pair<double,double>> boundaries;
//    // Bunny - geometry[1]
//    pair<double,double> x_limits(0., 1.);
//    pair<double,double> y_limits(0., 0.9);
//    pair<double,double> z_limits(0., 1.);
//    double gridspacing = 0.018*2;
//    double radius = 0.95*gridspacing;

    // funny - geometry[2]
    pair<double,double> x_limits(0., 1.);
    pair<double,double> y_limits(0., 0.9);
    pair<double,double> z_limits(0., 1.);
    double gridspacing = 0.008*2.5;
    double radius = 0.95*gridspacing;

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

    updateVerletList(verletStringId, particles, grid, radius);


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
        int gridId = grid.gridId(R.col(idCol.second));
        data(idCol.second, gridID) = gridId;
        data(idCol.second, connectionsID) = particles.verletList(idCol.first).size();
    }

    SaveParticles save_lmp("lmp");
    save_lmp.writeToFile(particles, string(TEST_SAVE_PATH) + "/verlet.lmp");
    SaveParticles save_xyz("xyz");
    save_xyz.writeToFile(particles, string(TEST_SAVE_PATH) + "/verlet.xyz");
}


