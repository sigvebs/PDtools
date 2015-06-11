#include <gtest/gtest.h>
#include <PdFunctions/pdfunctions.h>
#include <PDtools.h>
#include <PDtools/Solver/solvers.h>
#include <PDtools/Force/forces.h>

#include <PDtools/Modfiers/modifiers.h>

extern std::vector<string> geometries;

using namespace PDtools;

class PD_SOLVER_FIXTURE : public ::testing::Test {
    protected:

    PD_SOLVER_FIXTURE()
    {

    }
};

TEST_F(PD_SOLVER_FIXTURE, TEST_INTEGRATOR)
{
    using namespace std;
    using namespace arma;

//    double dt = 5000.65212887624e-10;
    double dt = 1.65212887624e-10/30.;
    int nSteps = 150;
    nSteps     = 10;

    int dim = 3;
    int saveFreq = 1500;
    double vMag = -0.5;
//    double s0 = 0.0158377370467;
    double s0 = 0.00158377370467;
    double rho = 8140.0;
    double alpha = 0.25;
    double delta = 3.66216165e-06;
    double gridspacing = 1.1*delta;

    double h = 1.25e-05;
    double E =  15200000000.0;
    double nu = 1./3.;
    double k = E/(2.*(1. - nu));
    double micromodulus = 12.*k/(M_PI*h*pow(delta, 3));
    micromodulus = k;

    vector<string> saveParameters = {"id", "x", "y", "stress", "damage"};
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------
    PD_Particles particles;
    particles = load_pd(geometries[8]);

    particles.registerParameter("rho", rho);
    particles.registerParameter("micromodulus", micromodulus);
    particles.registerParameter("s0", s0);
    particles.registerParameter("radius", s0);

    vector<pair<double,double>> boundaries;
    // Hole - geometry[4]
    //        pair<double,double> x_limits(0., 0.0487619);
    //        pair<double,double> y_limits(0., 0.0488042);
    //        pair<double,double> z_limits(0., 0.001);
    //        double gridspacing = 0.0022*3;;
    //        double radius = 0.95*gridspacing;

    //    // Bunny - geometry[1]
    //    pair<double,double> x_limits(0., 1.);
    //    pair<double,double> y_limits(0., 0.9);
    //    pair<double,double> z_limits(0., 1.);
    //    double gridspacing = 0.018*2;
    //    double radius = 0.95*gridspacing;

    // Bunny - geometry[1]
    pair<double,double> x_limits(0., 5e-05);
    pair<double,double> y_limits(0., 5e-05);
    pair<double,double> z_limits(-0.5*h, 0.5*h);

    boundaries = {x_limits, y_limits, z_limits};

    Grid grid(boundaries, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(particles);

    setPdConnections(particles, grid, delta);
    reCalculatePdMicromodulus(particles, dim);
    calculateRadius(particles, 2, h);

    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    Force * pdForce = new PD_bondForce(particles);
    Force * contactForce = new ContactForce(particles, grid, delta/3.01);
    Modifier *fracture= new PmbFracture(alpha);
    fracture->setParticles(particles);
    //--------------------------------------------------------------------------
    // Setting the boundary conditions
    //--------------------------------------------------------------------------
    double d = 0.4;
    pair<double, double> bound1(x_limits.first - d*delta, x_limits.first + d*delta);
    Modifier *boundaryLeft= new VelocityBoundary(-vMag, 0, bound1, 0, 500);
    boundaryLeft->setParticles(particles);

    pair<double, double> bound2(x_limits.second - d*delta, x_limits.second + d*delta);
    Modifier *boundaryRight = new VelocityBoundary(vMag, 0, bound2, 0, 500);
    boundaryRight->setParticles(particles);

    boundaryLeft->initialize();
    boundaryRight->initialize();

    double d_tb = 0.2;
    pair<double, double> boundTop(y_limits.first - d_tb*delta, y_limits.first + d_tb*delta);
    Modifier *boundaryTop = new VelocityBoundary(0, 1, boundTop, 1, 1);
    boundaryTop->setParticles(particles);

    pair<double, double> boundBottom(y_limits.second - d_tb*delta, y_limits.second + d_tb*delta);
    Modifier *boundaryBottom = new VelocityBoundary(0, 1, boundBottom, 1, 1);
    boundaryBottom->setParticles(particles);

    boundaryTop->initialize();
    boundaryBottom->initialize();

    fracture->initialize();
    //--------------------------------------------------------------------------
    // Setting the integrator
    //--------------------------------------------------------------------------
    Solver *solver = new VelocityVerletIntegrator();
    solver->setMainGrid(grid);
    solver->setParticles(particles);
    solver->addForce(pdForce);
    solver->addForce(contactForce);
    solver->setDt(dt);
    solver->setSteps(nSteps);
    solver->setSaveInterval(saveFreq);
    solver->setSaveParameters(saveParameters);
    solver->addModifier(boundaryLeft);
    solver->addModifier(boundaryRight);
    solver->addModifier(boundaryTop);
    solver->addModifier(boundaryBottom);
    solver->addSpModifier(fracture);
    solver->solve();

    delete solver;
    delete pdForce;
    delete boundaryLeft;
    delete boundaryRight;
    delete fracture;
}


TEST_F(PD_SOLVER_FIXTURE, TEST_ADR)
{
    using namespace std;
    using namespace arma;

    double dt = 1.0;
    int nSteps = 500;
    int dim = 2;

    int saveFreq = 100;
    double vMag = 1e-11;
    double s0 = 0.00158377370467;
    double rho = 2140.0;
    double alpha = 0.25;
    double delta = 3.66216165e-06;
    double gridspacing = 1.1*delta;

    double h = 1.25e-05;
    double E =  15200000000.0;
    double nu = 1./3.;
    double k = E/(2.*(1. - nu));
    double micromodulus = 12.*k/(M_PI*h*pow(delta, 3));
    micromodulus = k;

    double G0 = 0.1; // Fracture energy release rate

    vector<string> saveParameters = {"id", "x", "y", "potential_energy", "kinetic_energy", "stress", "damage"};
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------
    PD_Particles particles;
    particles = load_pd(geometries[8]);
    particles.initializeADR();
    particles.registerParameter("rho", rho);
    particles.registerParameter("micromodulus", micromodulus);
    particles.registerParameter("s0", s0);
    particles.registerParameter("radius", s0);

    vector<pair<double,double>> boundaries;

    // Bunny - geometry[1]
    pair<double,double> x_limits(0., 5e-05);
    pair<double,double> y_limits(0., 5e-05);
    pair<double,double> z_limits(-0.5*h, 0.5*h);

    boundaries = {x_limits, y_limits, z_limits};

    Grid grid(boundaries, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(particles);

    setPdConnections(particles, grid, delta);
    reCalculatePdMicromodulus(particles, dim);
    reCalculatePdFractureCriterion(particles, G0, delta, h);
    calculateRadius(particles, 2, h);
    //--------------------------------------------------------------------------
    // Setting the initial position
    //--------------------------------------------------------------------------
    mat &r0 = particles.r0();
    mat &r = particles.r();

    for(int i=0; i<particles.nParticles(); i++)
    {
        for(int d=0; d<3; d++)
        {
            r0(d, i) = r(d,i);
        }
    }
    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    Force * pdForce = new PD_bondForce(particles);
    Force * contactForce = new ContactForce(particles, grid, delta/2.01);
    Modifier *fracture= new PmbFracture(alpha);
    fracture->setParticles(particles);
    //--------------------------------------------------------------------------
    // Setting the boundary conditions
    //--------------------------------------------------------------------------
    double d = 0.75;
    pair<double, double> bound1(x_limits.first - d*delta, x_limits.first + d*delta);
    Modifier *boundaryLeft= new MoveParticles(-vMag, 0, bound1, 0);
    boundaryLeft->setParticles(particles);

    pair<double, double> bound2(x_limits.second - d*delta, x_limits.second + d*delta);
    Modifier *boundaryRight = new MoveParticles(vMag, 0, bound2, 0);
    boundaryRight->setParticles(particles);

    boundaryLeft->initialize();
    boundaryRight->initialize();

    fracture->initialize();

    //--------------------------------------------------------------------------
    // Testing the surface scaling
    //--------------------------------------------------------------------------
    vector<Force*> forces = {pdForce};
    surfaceCorrection(particles, forces, k, nu, dim);

    //--------------------------------------------------------------------------
    // Setting the integrator
    //--------------------------------------------------------------------------
    Solver *solver = new ADR();
    solver->setMainGrid(grid);
    solver->setParticles(particles);
    solver->addForce(pdForce);
    solver->addForce(contactForce);
    solver->setDt(dt);
    solver->setSteps(nSteps);
    solver->setSaveInterval(saveFreq);
    solver->setSaveParameters(saveParameters);
    solver->addModifier(boundaryLeft);
    solver->addModifier(boundaryRight);
    solver->addSpModifier(fracture);
    solver->solve();

    delete solver;
    delete pdForce;
    delete boundaryLeft;
    delete boundaryRight;
}


TEST_F(PD_SOLVER_FIXTURE, TEST_ADR_BUNNY)
{
    using namespace std;
    using namespace arma;

    double dt = 1.0;
    int nSteps = 10;
    int dim = 3;

    int saveFreq = 10;
    double vMag = -1e-11;
    double s0 = 0.00158377370467;
    double rho = 2140.0;

    double alpha = 0.25;
    double delta = 0.018*2;
    double gridspacing = 1.1*delta;

    double E =  15200000000.0;
    double nu = 1./3.;
    double k = E/(2.*(1. - nu));
    double micromodulus = k;
    double error = 1e-4;

    double G0 = 0.1; // Fracture energy release rate

    vector<string> saveParameters = {"id", "x", "y", "z", "potential_energy", "kinetic_energy", "stress", "damage"};
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------
    PD_Particles particles;
    particles = load_pd(geometries[1]);
    particles.initializeADR();
    particles.registerParameter("rho", rho);
    particles.registerParameter("micromodulus", micromodulus);
    particles.registerParameter("s0", s0);
    particles.registerParameter("radius", s0);

    vector<pair<double,double>> boundaries;

    // Bunny - geometry[1]
    pair<double,double> x_limits(0., 1.);
    pair<double,double> y_limits(0., 0.9);
    pair<double,double> z_limits(0., 1.);

    boundaries = {x_limits, y_limits, z_limits};

    Grid grid(boundaries, gridspacing);
    grid.initialize();
    grid.placeParticlesInGrid(particles);

    setPdConnections(particles, grid, delta);
    reCalculatePdMicromodulus(particles, dim);
    reCalculatePdFractureCriterion(particles, G0, delta);
    calculateRadius(particles, 3);
    //--------------------------------------------------------------------------
    // Setting the initial position
    //--------------------------------------------------------------------------
    mat &r0 = particles.r0();
    mat &r = particles.r();

    for(int i=0; i<particles.nParticles(); i++)
    {
        for(int d=0; d<3; d++)
        {
            r0(d, i) = r(d,i);
        }
    }
    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    Force * pdForce = new PD_bondForce(particles);
    Force * contactForce = new ContactForce(particles, grid, delta/1.21);
    Modifier *fracture= new PmbFracture(alpha);
    fracture->setParticles(particles);
    //--------------------------------------------------------------------------
    // Setting the boundary conditions
    //--------------------------------------------------------------------------
    double d = 1.3;

    pair<double, double> bound1(x_limits.first - d*delta, x_limits.first + d*delta);
    Modifier *boundaryLeft= new MoveParticles(-vMag, 0, bound1, 0);
    boundaryLeft->setParticles(particles);

    pair<double, double> bound2(x_limits.second - 2*d*delta, x_limits.second + d*delta);
    Modifier *boundaryRight = new MoveParticles(vMag, 0, bound2, 0);
    boundaryRight->setParticles(particles);

    boundaryLeft->initialize();
    boundaryRight->initialize();

    fracture->initialize();
    //--------------------------------------------------------------------------
    // Setting the integrator
    //--------------------------------------------------------------------------
    Solver *solver = new ADR();
    solver->setErrorThreshold(error);
    solver->setMainGrid(grid);
    solver->setParticles(particles);
    solver->addForce(pdForce);
    solver->addForce(contactForce);
    solver->setDt(dt);
    solver->setSteps(nSteps);
    solver->setSaveInterval(saveFreq);
    solver->setSaveParameters(saveParameters);
    solver->addModifier(boundaryLeft);
    solver->addModifier(boundaryRight);
    solver->addSpModifier(fracture);
    solver->solve();

    delete solver;
    delete pdForce;
    delete boundaryLeft;
    delete boundaryRight;
}
