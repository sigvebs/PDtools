#include <gtest/gtest.h>
#include <memory>
#include <libconfig.h++>

#include <PDtools.h>

using namespace libconfig;
using namespace PDtools;

extern std::vector<string> configurations;

class CONFIGURATION_FIXTURE : public ::testing::Test {
    protected:
    PDtools::Particles testParticles;

    CONFIGURATION_FIXTURE()
    {

    }
};

TEST_F(CONFIGURATION_FIXTURE, LOAD_PD_CONFIGURATION)
{
    using namespace std;
    Config cfg;
    string cfgFileName = configurations[0];
    cout << cfgFileName << endl;
    cfg.readFile(cfgFileName.c_str());
    Setting & root = cfg.getRoot();

    PDsharedData PDdata;
    //--------------------------------------------------------------------------
    // Loading the domain
    //--------------------------------------------------------------------------
    const Setting & cfg_domain = root["domain"];
    const Setting & cfg_boundaries = cfg_domain["boundaries"];
    const Setting & cfg_periodicity = cfg_domain["periodic"];

    pair<double,double> x_limits(cfg_boundaries[0], cfg_boundaries[1]);
    pair<double,double> y_limits(cfg_boundaries[2], cfg_boundaries[3]);
    pair<double,double> z_limits(cfg_boundaries[4], cfg_boundaries[5]);
    vector<std::pair<double,double>> boundaries = {x_limits, y_limits, z_limits};
    ivec3 periodicBoundary = {cfg_periodicity[0], cfg_periodicity[1], cfg_periodicity[2]};
    int dim = cfg_domain["dim"];

    Domain *domain = new Domain(dim, boundaries);
    domain->periodicBoundaries(periodicBoundary);
    PDdata.domain(domain);

    //--------------------------------------------------------------------------
    // Setting the Peridynamic grid
    //--------------------------------------------------------------------------
    double gridSpacing = cfg_domain["spacing"];
    Grid *pdGrid = new Grid(PDdata.domain(), gridSpacing);
    pdGrid->update();
    pdGrid->createGrid();
    PDdata.pdGrid(pdGrid);

    //--------------------------------------------------------------------------
    // Loading the particles from file
    //--------------------------------------------------------------------------
    const Setting & cfg_particles = root["particles"];
    for(int i=0; i<cfg_particles.getLength(); i++)
    {
        const Setting & particleSetting = cfg_particles[i];
//        string type = particleSetting["type"];
        string particlePath = particleSetting["particlePath"];
        unordered_map<string, double> particleParameters;

        const Setting & pm = particleSetting["particleParameters"];
        for(int j=0;j<pm.getLength(); j++)
        {
            particleParameters[pm[j].getName()] = pm[j];
        }

        PDdata.particles(new PD_Particles(load_pd(particlePath, particleParameters)));
    }

    //--------------------------------------------------------------------------
    // Setting the integrator
    //--------------------------------------------------------------------------
    PDdata.t(0);
    PDdata.timeStep(1.0);
    PDdata.nIterations(20);
    PDdata.iteration(0);
    //--------------------------------------------------------------------------
}
