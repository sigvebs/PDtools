#include <gtest/gtest.h>
#include <libconfig.h++>
#include <PdFunctions/pdfunctions.h>
#include <PDtools.h>
#include <PDtools/Solver/solvers.h>
#include <PDtools/Force/forces.h>
#include "PDtools/SavePdData/savepddata.h"

#include <PDtools/Modfiers/modifiers.h>

extern std::vector<string> geometries;

using namespace PDtools;

class PD_STRETCH_FIXTURE : public ::testing::Test {
    protected:

    PD_STRETCH_FIXTURE()
    {

    }
};


void test_stretch(string project_folder, int n, string save_folder)
{
    double epsilon = 0.05;
    using namespace std;
    using namespace arma;
    using namespace libconfig;

    string cfgPath = project_folder + "/" + to_string(n) + "/run.cfg";
    cout << cfgPath << endl;
    Config m_cfg;

    try
    {
        m_cfg.readFile(cfgPath.c_str());
    }
    catch (const FileIOException &fioex)
    {
        std::cerr << "I/O error while reading the configuration file." << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        exit(EXIT_FAILURE);
    }
    m_cfg.setAutoConvert(true);

    int dim;
    double s0;
    double rho;
    double delta;
    double h;
    double E;
    double nu;
    double k;

    m_cfg.lookupValue("dim", dim);
    m_cfg.lookupValue("E", E);
    m_cfg.lookupValue("nu", nu);
    m_cfg.lookupValue("delta", delta);
    m_cfg.lookupValue("rho", rho);
    double gridspacing = 1.1*delta;

    string particles_path;
    m_cfg.lookupValue("particlesPath", particles_path);

    vector<string> saveParameters = {"id", "x", "y", "potential_energy", "kinetic_energy", "stress", "max_stretch"};

    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    PD_Particles particles;
    particles = load_pd(particles_path);

    particles.initializeADR();
//    particles.registerParameter("s0", s0);

    Setting &cfg_domain = m_cfg.lookup("domain");

    double dxdydz[3];
    vector<pair<double,double>> domain;
    for(int d=0; d<3; d++)
    {
        pair<double, double> bound(cfg_domain[2*d], cfg_domain[2*d+1]);
        dxdydz[d] = bound.second - bound.first;
        domain.push_back(bound);
    }

    Grid grid(domain, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(particles);

    setPdConnections(particles, grid, delta, gridspacing);

    //--------------------------------------------------------------------------
    // Setting the initial position
    //--------------------------------------------------------------------------
    mat &r0 = particles.r0();
    mat &r = particles.r();

    for(int i=0; i<particles.nParticles(); i++)
    {
        for(int d=0; d<3; d++)
        {
            r0(i, d) = r(i, d);
        }
    }

    //--------------------------------------------------------------------------
    // Setting particles parameters
    //--------------------------------------------------------------------------

    h = dxdydz[2];
    double c;

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
    else
    {
        cerr << "ERROR: dimension " << dim << " not supported" << endl;
        cerr << "use 2 or 3." << endl;
        exit(EXIT_FAILURE);
    }
    bool useS0fromCfg = false;
    if(m_cfg.lookupValue("s0", s0))
    {
        useS0fromCfg = true;
    }

    particles.registerParameter("rho", rho);
    particles.registerParameter("radius");

    bool calculateMicromodulus = m_cfg.lookupValue("calculateMicromodulus", calculateMicromodulus);

    if(calculateMicromodulus){
        particles.registerParameter("micromodulus", k);
        reCalculatePdMicromodulus(particles, dim);
    }
    else
    {
        particles.registerParameter("micromodulus", c);
    }
    calculateRadius(particles, dim, dxdydz[2]);
    reCalculatePdMicromodulus(particles, dim);

    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    Force * pdForce = new PD_bondForce(particles);

    //--------------------------------------------------------------------------
    // Testing the surface scaling
    //--------------------------------------------------------------------------
    vector<Force*> forces = {pdForce};
    surfaceCorrection(particles, forces, k, nu, dim);

    //--------------------------------------------------------------------------
    // Updating the stretch
    //--------------------------------------------------------------------------
    string save_path = save_folder;
    SavePdData *saveParticles = new SavePdData(saveParameters);
    saveParticles->setSavePath(save_path);
    saveParticles->setParticles(&particles);
    saveParticles->setForces(forces);
    saveParticles->initialize();

    //--------------------------------------------------------------------------
    // Stretching the system
    //--------------------------------------------------------------------------
    for(int i=0; i<particles.nParticles(); i++)
    {
        r(i, 0) *= (1 + epsilon);
        r(i, 1) *= (1 - nu*epsilon);
    }
    //--------------------------------------------------------------------------
    // Calculating the forces
    //--------------------------------------------------------------------------
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<particles.nParticles(); i++)
    {
        pair<int, int> id(i, i);
        pdForce->calculateForces(id);
    }

    saveParticles->evaluate(1, 1);
    //--------------------------------------------------------------------------
    // Saving the results
    //--------------------------------------------------------------------------
    SaveParticles save_xyz("lmp");
    save_path = save_folder + "/" + to_string(n) + ".lmp";
    save_xyz.writeToFile(particles, save_path);
}

TEST_F(PD_STRETCH_FIXTURE, TEST_STRETCH)
{
    string project_folder = "/media/sigve/Pengebingen/scratch/Projects/PlateWithHole/adrConstantConvergence";
    string save_folder = "/media/sigve/Pengebingen/scratch/tmp/analysis";
    vector<int> n_particles = {2000, 4000, 6000, 8000, 10000, 20000, 40000, 60000, 80000, 100000};
    for(int n:n_particles)
    {
        test_stretch(project_folder, n, save_folder);
    }
}
