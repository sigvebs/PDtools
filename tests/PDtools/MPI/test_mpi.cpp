#include <gtest/gtest.h>

#include <PDtools.h>

using namespace PDtools;
using namespace arma;

extern std::vector<string> geometries;
//------------------------------------------------------------------------------
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

        r = _r.t();
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
//------------------------------------------------------------------------------
