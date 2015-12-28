#include "pdsolver.h"

#include <PDtools/PdFunctions/pdfunctions.h>
#include <PDtools/Force/forces.h>
#include <PDtools/Modfiers/modifiers.h>
#include <PDtools/Solver/solvers.h>

#include <armadillo>
#include <boost/regex.hpp>
#include <unordered_map>
#include <vector>


using namespace PDtools;
using namespace arma;
using namespace libconfig;
//------------------------------------------------------------------------------
PdSolver::PdSolver(string cfgPath):
    m_configPath(cfgPath)
{
    // Reading configuration file
    try
    {
        m_cfg.readFile(cfgPath.c_str());
    }
    catch (const FileIOException &fioex)
    {
        std::cerr << "I/O error while reading the configuration file."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const ParseException &pex)
    {
        std::cerr << "Parse error at "
                  << pex.getFile()
                  << ":"
                  << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        exit(EXIT_FAILURE);
    }
    m_cfg.setAutoConvert(true);

    // TODO: TMP FIXES, WEIRD COMPILER STUFF
//    vec be(2);
//    be.randu();
//    boost::regex rr("([A-Za-z_]+)");
}
//------------------------------------------------------------------------------
PdSolver::~PdSolver()
{
}
//------------------------------------------------------------------------------
void PdSolver::initialize()
{
    int dim = 3;
    double E = -1;
    double nu = -1;
    double G0 = -1;
    double rho = -1;

    double E0 = 1;
    double L0 = 1;
    double v0 = 1;
    double t0 = 1;
    double rho0 = 1;

    m_cfg.lookupValue("dim", dim);
    m_cfg.lookupValue("E", E);
    m_cfg.lookupValue("nu", nu);
    m_cfg.lookupValue("G0", G0);
    m_cfg.lookupValue("rho", rho);

    // Checking for dimensional scaling
    int scaleparameters = 0;
    m_cfg.lookupValue("scaleparameters", scaleparameters);

    if(scaleparameters)
    {
        m_cfg.lookupValue("E0", E0);
        m_cfg.lookupValue("L0", L0);
        m_cfg.lookupValue("v0", v0);
        m_cfg.lookupValue("t0", t0);
        m_cfg.lookupValue("rho0", rho0);
    }

    //--------------------------------------------------------------------------
    // Setting the domain
    //--------------------------------------------------------------------------
    if(!m_cfg.exists("domain"))
    {
        cerr << "'domain' must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    Setting &cfg_domain = m_cfg.lookup("domain");

    if(cfg_domain.getLength() != 6)
    {
        cerr << "'domain' must be set on the form ";
        cerr << "'domain = [x0, x1, y0, y1, z0, z1]'" << endl;
        exit(EXIT_FAILURE);
    }

    double dxdydz[3];
    vector<pair<double,double>> domain;
    for(int d=0; d<3; d++)
    {
        pair<double, double> bound(cfg_domain[2*d], cfg_domain[2*d+1]);
        bound.first /= L0;
        bound.second /= L0;
        dxdydz[d] = bound.second - bound.first;
        domain.push_back(bound);
    }
    // TODO: periodic boundaries
    //--------------------------------------------------------------------------
    // Loading the particles
    //--------------------------------------------------------------------------
    if(!m_cfg.exists("particlesPath"))
    {
        cerr << "'particlesPath' must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    string particlesPath = (const char *) m_cfg.lookup("particlesPath");
    m_particles = load_pd(particlesPath);
    m_particles.dimensionalScaling(E0, L0, v0, t0, rho0);
    cout << "Loading particles complete" << endl;
    //--------------------------------------------------------------------------
    // TODO: Setting the initial position, should not be done here
    //--------------------------------------------------------------------------
    mat &r0 = m_particles.r0();
    const mat &r = m_particles.r();

    for(unsigned int i=0; i<m_particles.nParticles(); i++)
    {
        for(int d=0; d<3; d++)
        {
            r0(i, d) = r(i, d);
        }
    }
    //--------------------------------------------------------------------------
    // Setting the grid
    //--------------------------------------------------------------------------
    double delta;
    if(!m_cfg.lookupValue("delta", delta))
    {
        cerr << "'delta' must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    double lc;
    if(!m_cfg.lookupValue("lc", lc))
    {
        cerr << "'lc' must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    lc /= L0;
    delta /= L0;

    cout << "gridspacing" << endl;
    const double gridspacing = 1.2*delta;
    cout << "Grid" << endl;
    m_grid = Grid(domain, gridspacing);
    m_grid.dim = dim;
    cout << "initialize" << endl;
    m_grid.initialize();
    cout << "placeParticlesInGrid" << endl;
    m_grid.placeParticlesInGrid(m_particles);
    cout << "Particles placed in grid" << endl;
    //--------------------------------------------------------------------------
    // Setting particles parameters
    //--------------------------------------------------------------------------
    double h = dxdydz[2];
    double s0 = 1;

    bool useS0fromCfg = false;
    if(m_cfg.lookupValue("s0", s0))
    {
        useS0fromCfg = true;
    }

    m_particles.registerParameter("rho", rho/rho0);
    m_particles.registerParameter("s0", s0);
    m_particles.registerParameter("radius");
    calculateRadius(m_particles, dim, dxdydz[2]);

    int calculateMicromodulus = 0;
    m_cfg.lookupValue("calculateMicromodulus", calculateMicromodulus);
    cout << "Additional parameters set" << endl;
    //--------------------------------------------------------------------------
    // Setting the PD-connections
    //--------------------------------------------------------------------------
    //lc *= 1.05;

    bool applyVolumeCorrectio = true;
    setPdConnections(m_particles, m_grid, delta, lc);
    addFractures(m_particles, domain);
    removeVoidConnections(m_particles, m_grid, delta, lc);
    cleanUpPdConnections(m_particles);
    m_particles.registerPdParameter("volumeScaling", 1);

    if(applyVolumeCorrectio)
    {
        applyVolumeCorrection(m_particles, delta, lc);
    }
    setPD_N3L(m_particles);
    cout << "PD connections set" << endl;
    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    vector<Force*> forces;

    Setting &cfg_forces = m_cfg.lookup("forces");
    for(int i=0; i<cfg_forces.getLength(); i++)
    {
        const char * tmpType;
        cfg_forces[i].lookupValue("type", tmpType);
        string type = tmpType;
        delete tmpType;

        if(type == "bond force")
        {
            forces.push_back(new PD_bondForce(m_particles));
        }
        else if(type == "LPS")
        {
            forces.push_back(new PD_LPS(m_particles));
        }
        else if(type == "OSP")
        {
            forces.push_back(new PD_OSP(m_particles));
        }
        else if(type == "PMB")
        {
            double alpha;

            if (!cfg_forces[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for force '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            forces.push_back(new PD_PMB(m_particles, lc, delta, alpha));
        }
        else if(type == "gaussian bond force")
        {
            forces.push_back(new PD_bondforceGaussian(m_particles));
        }
        else if(type == "weighted bond")
        {
            string weightType;
            cfg_forces[i].lookupValue("weightType", weightType);
            forces.push_back(new PD_bondforceGaussian(m_particles, weightType));
        }
        else if(type == "contact force")
        {
            forces.push_back(new ContactForce(m_particles, m_grid, lc));
        }
    }


    // Initializing forces
    for(Force* force:forces)
    {
        force->numericalInitialization(calculateMicromodulus);
        force->initialize(E/E0, nu, delta, dim, h, lc);
    }
    cout << "Forces set" << endl;
    //--------------------------------------------------------------------------
    // Recalcuating particle properties
    //--------------------------------------------------------------------------
    int applySurfaceCorrection = 0;
    m_cfg.lookupValue("applySurfaceCorrection", applySurfaceCorrection);
    if(applySurfaceCorrection)
    {
        for(Force* force:forces)
        {
            int nSurfaceCorrections = 20;
            for(int i=0; i<nSurfaceCorrections; i++)
            {
                //force->applySurfaceCorrection(0.0005);
                force->applySurfaceCorrection(0.001);
            }
        }
    }

    if(!useS0fromCfg)
    {
        if(dim == 3)
        {
            reCalculatePdFractureCriterion(m_particles, G0, delta);
        }
        if(dim == 2)
        {
            reCalculatePdFractureCriterion(m_particles, G0, delta, h);
        }
        if(dim == 1)
        {
            cerr << "ERROR: dimension " << dim << " not supported for numerical 's0'" << endl;
            cerr << "use 2 or 3." << endl;
            exit(EXIT_FAILURE);
        }
    }
    cout << "Corrections set" << endl;
    //--------------------------------------------------------------------------
    // Setting the solver
    //--------------------------------------------------------------------------
    string solverType = (const char *) m_cfg.lookup("solverType");

    int nSteps;
    double dt;

    if (!m_cfg.lookupValue("nSteps", nSteps))
    {
        cerr << "Error reading the 'nSteps' in config file" << endl;
        exit(EXIT_FAILURE);
    }

    if(solverType == "ADR")
    {
        m_particles.initializeADR();
        ADR *adrSolver = new ADR();
        dt = 1.0;
        double errorThreshold;
        int maxSteps = 2000;
        int maxStepsFracture = 1000;
        m_cfg.lookupValue("maxSteps", maxSteps);
        m_cfg.lookupValue("maxStepsFracture", maxStepsFracture);

        if(!m_cfg.lookupValue("errorThreshold", errorThreshold))
        {
            cerr << "Error reading the 'errorThreshold' in config file" << endl;
            exit(EXIT_FAILURE);
        }

        adrSolver->maxSteps(maxSteps);
        adrSolver->maxStepsFracture(maxStepsFracture);
        adrSolver->setErrorThreshold(errorThreshold);
        solver = adrSolver;
    }
    else if(solverType == "dynamic ADR")
    {
        m_particles.initializeADR();
        solver = new dynamicADR();
        m_cfg.lookupValue("dt", dt);
    }
    else if(solverType == "conjugate gradient")
    {
        const int maxIterations = 20000;
        const int threshold = 3.e-8;
        solver = new StaticSolver(maxIterations, threshold);
        m_cfg.lookupValue("dt", dt);
    }
    else if(solverType == "velocity verlet")
    {
        solver = new VelocityVerletIntegrator();
        m_cfg.lookupValue("dt", dt);
    }
    else if(solverType == "euler-chromer")
    {
        solver = new EulerCromerIntegrator();
        m_cfg.lookupValue("dt", dt);
    }
    else
    {
        cerr << "Error: solver not set" << endl;
        exit(EXIT_FAILURE);
    }
    dt /= t0;

    cout << "Solver set" << endl;
    //--------------------------------------------------------------------------
    // Setting the modifiers
    //--------------------------------------------------------------------------
    vector<Modifier *> modifiers;
    vector<Modifier *> spModifiers;
    vector<Modifier *> qsModifiers;

    Setting &cfg_modifiers = m_cfg.lookup("modifiers");
    for(int i=0; i<cfg_modifiers.getLength(); i++)
    {
        const char * tmpType;
        cfg_modifiers[i].lookupValue("type", tmpType);
        string type = tmpType;
        delete tmpType;

        if(type == "velocity particles")
        {
            double v;
            int axis, vAxis;
            int isStatic = 0;
            int nSteps = 1;

            double a0 = cfg_modifiers[i]["area"][0];
            double a1 = cfg_modifiers[i]["area"][1];
            cfg_modifiers[i].lookupValue("static", isStatic);
            a0 /= L0;
            a1 /= L0;

            pair<double, double> area(a0, a1);

            if (!cfg_modifiers[i].lookupValue("v", v) ||
                    !cfg_modifiers[i].lookupValue("axis", axis) ||
                    !cfg_modifiers[i].lookupValue("vAxis", vAxis))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            v /= v0;
            modifiers.push_back(new VelocityBoundary(v, vAxis, area, axis,
                                                     dt, nSteps, isStatic));
        }
        else if(type == "move particles")
        {
            double v;
            int axis, vAxis;
            int isStatic = 0;

            double a0 = cfg_modifiers[i]["area"][0];
            double a1 = cfg_modifiers[i]["area"][1];
            cfg_modifiers[i].lookupValue("static", isStatic);
            a0 /= L0;
            a1 /= L0;

            pair<double, double> area(a0, a1);

            if (!cfg_modifiers[i].lookupValue("v", v) ||
                    !cfg_modifiers[i].lookupValue("axis", axis) ||
                    !cfg_modifiers[i].lookupValue("vAxis", vAxis))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            v /= v0;

            modifiers.push_back(new MoveParticles(v, vAxis, area, axis, dt, isStatic));
        }
        else if(type == "boundary force")
        {
            double appliedForce;
            int axis, forceAxis;

            double a0 = cfg_modifiers[i]["area"][0];
            double a1 = cfg_modifiers[i]["area"][1];
            a0 /= L0;
            a1 /= L0;


            pair<double, double> area(a0, a1);

            if (!cfg_modifiers[i].lookupValue("appliedForce", appliedForce) ||
                    !cfg_modifiers[i].lookupValue("axis", axis) ||
                    !cfg_modifiers[i].lookupValue("forceAxis", forceAxis))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            appliedForce /= (E0/pow(L0, 4));

            if(solverType == "ADR")
            {
                modifiers.push_back(new boundaryForce(appliedForce, forceAxis, area, axis));
            }
            else
            {
                modifiers.push_back(new boundaryForce(appliedForce, forceAxis, area, axis));
            }
        }
        else if(type == "PMB fracture")
        {
            double alpha;

            if (!cfg_modifiers[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            if(solverType == "ADR")
            {
                qsModifiers.push_back(new ADRfracture(alpha));
            }
            else
            {
                spModifiers.push_back(new PmbFracture(alpha));
            }
        }
        else if(type == "bond energy fracture")
        {
            double G;

            if (!cfg_modifiers[i].lookupValue("G", G))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            Modifier * fMod = new BondEnergyFracture(delta, G, dim, forces, h);
            if(solverType == "ADR")
            {
                qsModifiers.push_back(fMod);
            }
            else
            {
                spModifiers.push_back(fMod);
            }
        }
        else if(type == "Mohr-Coulomb fracture")
        {
            double mu, C, T;

            if (!cfg_modifiers[i].lookupValue("mu", mu)||
                    !cfg_modifiers[i].lookupValue("C", C) ||
                    !cfg_modifiers[i].lookupValue("T", T) )
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            C /= E0;
            T /= E0;

            if(solverType == "ADR")
            {
                ADRmohrCoulombFracture * failureCriterion = new ADRmohrCoulombFracture(mu, C, T, dim);
                for(Force* force: forces)
                {
                    failureCriterion->addForce(force);
                }
                qsModifiers.push_back(failureCriterion);
            }
            else
            {
                MohrCoulombFracture * failureCriterion = new MohrCoulombFracture(mu, C, T, dim);
                for(Force* force: forces)
                {
                    failureCriterion->addForce(force);
                }
                spModifiers.push_back(failureCriterion);
            }
        }
        else if(type == "simple fracture")
        {
            double alpha;

            if (!cfg_modifiers[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            SimpleFracture * simpleFracture = new SimpleFracture(alpha);
            if(solverType == "ADR")
            {
                qsModifiers.push_back(simpleFracture);
            }
            else
            {
                spModifiers.push_back(simpleFracture);
            }
        }
        else if(type == "ADR fracture average")
        {
            double alpha;

            if (!cfg_modifiers[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            ADRfractureAverage * adrFracture = new ADRfractureAverage(alpha);
            qsModifiers.push_back(adrFracture);
        }
        else
        {
            cerr << "ERROR: modifier '" << type << "' does not exsits." << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Initializing modifiers
    for(Modifier* mod: modifiers)
    {
        mod->setParticles(m_particles);
        mod->initialize();
    }

    for(Modifier* mod: spModifiers)
    {
        mod->setParticles(m_particles);
        mod->initialize();
    }

    for(Modifier* mod: qsModifiers)
    {
        mod->setParticles(m_particles);
        mod->initialize();
    }

    cout << "Modifiers set" << endl;
    //--------------------------------------------------------------------------
    // Setting additional initial conditions
    //--------------------------------------------------------------------------
    if( m_cfg.exists("initialConditions"))
    {
        Setting &cfg_initialConditions = m_cfg.lookup("initialConditions");
        for(int i=0; i<cfg_initialConditions.getLength(); i++)
        {
            const char * tmpType;
            cfg_initialConditions[i].lookupValue("type", tmpType);
            string type = tmpType;
            delete tmpType;

            if(type == "strain")
            {
                double strain;
                int axis;

                double a0 = cfg_initialConditions[i]["area"][0];
                double a1 = cfg_initialConditions[i]["area"][1];
                a0 /= L0;
                a1 /= L0;

                pair<double, double> area(a0, a1);

                if (!cfg_initialConditions[i].lookupValue("strain", strain) ||
                        !cfg_initialConditions[i].lookupValue("axis", axis) ||
                        !cfg_initialConditions[i].lookupValue("axis", axis))
                {
                    cerr << "Error reading the parameters for modifier '"
                         << type << "'" << endl;
                    exit(EXIT_FAILURE);
                }
                applyInitialStrainStrain(m_particles, strain, axis, area);
            }
        }
    }

    cout << "Initial conditions set" << endl;
    //--------------------------------------------------------------------------
    // Adding the modifiers and forces to the solver
    //--------------------------------------------------------------------------
    solver->setMainGrid(m_grid);
    solver->setParticles(m_particles);
    solver->setSteps(nSteps);
    solver->setDt(dt);

    for(Force* force: forces)
    {
        solver->addForce(force);
    }
    for(Modifier* mod: modifiers)
    {
        solver->addModifier(mod);
    }
    for(Modifier* mod: spModifiers)
    {
        solver->addSpModifier(mod);
    }
    for(Modifier* mod: qsModifiers)
    {
        solver->addQsModifiers(mod);
    }

    //--------------------------------------------------------------------------
    // Setting the save values
    //--------------------------------------------------------------------------
    int saveFrequency;
    if (!m_cfg.lookupValue("saveFrequency", saveFrequency))
    {
        cerr << "Error reading the 'saveFrequency' in config file" << endl;
        exit(EXIT_FAILURE);
    }
    if(!m_cfg.exists("savePath"))
    {
        cerr << "Error reading the 'savePath' in config file" << endl;
        exit(EXIT_FAILURE);
    }
    string savePath = (const char *) m_cfg.lookup("savePath");

    vector<string> saveParameters;
    if(m_cfg.exists("saveParameters"))
    {
        const Setting & cfg_saveParameters = m_cfg.lookup("saveParameters");
        for(int i=0; i<cfg_saveParameters.getLength(); i++)
        {
            saveParameters.push_back(cfg_saveParameters[i].c_str());
        }
    }
    else
    {
        cerr << "No saveParamters set" << endl;
        exit(EXIT_FAILURE);
    }
    solver->setDim(dim);
    solver->setSavePath(savePath);
    solver->setSaveScaling(E0, L0, v0, t0, rho0);
    solver->setSaveInterval(saveFrequency);
    solver->setSaveParameters(saveParameters);
}
//------------------------------------------------------------------------------
void PdSolver::solve()
{
    cout << "Starting solver" << endl;
    solver->solve();
}
//------------------------------------------------------------------------------
