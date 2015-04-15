#include <gtest/gtest.h>
#include <PDtools.h>
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

        r = _r;
        v = _v;
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

        ASSERT_EQ(r_test(0, id_test), r(0, id_load));
        ASSERT_EQ(r_test(1, id_test), r(1, id_load));
        ASSERT_EQ(v_test(0, id_test), v(0, id_load));
        ASSERT_EQ(v_test(1, id_test), v(1, id_load));

        for(auto param:parameters)
        {
            int particles_type_pos = param.second;
            int testParticles_type_pos = parameters_test[param.first];
            ASSERT_EQ(data_test(id_test, testParticles_type_pos), data(id_load, particles_type_pos))
                    << "The data in Testparticles does not match the data read from file";
        }
    }
}
