#include "pdsolver.h"

#include <PDtools/PdFunctions/pdfunctions.h>
#include <PDtools/Force/forces.h>
#include <PDtools/Modfiers/modifiers.h>
#include <PDtools/Solver/solvers.h>
#include <PDtools/CalculateProperties/calculateproperties.h>
#include "PDtools/SavePdData/savepddata.h"

#include <armadillo>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace PDtools;
using namespace arma;
using namespace libconfig;
//------------------------------------------------------------------------------
PdSolver::PdSolver(string cfgPath, int myRank, int nMpiNodes):
    m_configPath(cfgPath),
    m_myRank(myRank),
    m_nCores(nMpiNodes)
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
    if(myRank == 0)
    {
        isRoot = true;
    }

}
//------------------------------------------------------------------------------
PdSolver::~PdSolver()
{
}
//------------------------------------------------------------------------------
int PdSolver::initialize()
{
    wall_clock timer;
    timer.tic();

    vector<pair<string, int>> neededProperties;

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
    if(isRoot)
        cout << "Setting the domain" << endl;

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
    // Setting the grid
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Setting the grid" << endl;
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

    double gridspacing = 1.15*(delta + 0.5*lc);
    m_grid = Grid(domain, gridspacing);
    m_grid.setIdAndCores(m_myRank, m_nCores);
    m_grid.dim(dim);
    m_grid.initialize();
    m_grid.setMyGridpoints();
    m_grid.setInitialPositionScaling(L0);
    //--------------------------------------------------------------------------
    // Loading the particles
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Loading the particles" << endl;
    if(!m_cfg.exists("particlesPath"))
    {
        cerr << "'particlesPath' must be set in the configuration file" << endl;
        exit(EXIT_FAILURE);
    }
    string particlesPath = (const char *) m_cfg.lookup("particlesPath");
    m_particles = load_pd(particlesPath, m_grid);
    m_particles.dimensionalScaling(E0, L0, v0, t0, rho0);
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
    // Setting particles parameters
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Setting particles parameters" << endl;
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
    //--------------------------------------------------------------------------
    // Setting the PD-connections
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Gridding particles and setting PD-connections" << endl;

    m_grid.clearParticles();
    m_grid.placeParticlesInGrid(m_particles);
    //lc *= 1.05;
#ifdef USE_MPI
    vector<string> ghostParameters = {"volume", "radius"};
    for(const string &param:ghostParameters)
    {
        m_particles.addGhostParameter(param);
    }
    exchangeGhostParticles(m_grid, m_particles);
#endif
    bool performVolumeCorrection = true;
    setPdConnections(m_particles, m_grid, delta, lc);
    m_grid.clearGhostParticles();
    exchangeInitialGhostParticles(m_grid, m_particles);
//    addFractures(m_particles, domain);
    removeVoidConnections(m_particles, m_grid, delta, lc);
    cleanUpPdConnections(m_particles);
    m_particles.registerPdParameter("volumeScaling", 1);
    if(performVolumeCorrection)
    {
        applyVolumeCorrection(m_particles, delta, lc);
    }

    setPD_N3L(m_particles);
    //--------------------------------------------------------------------------
    // Setting the Forces
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Setting the Forces" << endl;
    vector<string> forcesSet;
    vector<Force*> forces;
    Setting &cfg_forces = m_cfg.lookup("forces");

    for(int i=0; i<cfg_forces.getLength(); i++)
    {
        const char * tmpType;
        cfg_forces[i].lookupValue("type", tmpType);
        string type = tmpType;

        if(boost::iequals(type, "bond force"))
        {
            forces.push_back(new PD_bondForce(m_particles));
        }
        else if(boost::iequals(type, "LPS"))
        {
            forces.push_back(new PD_LPS(m_particles));
        }
        else if(boost::iequals(type, "OSP"))
        {
            forces.push_back(new PD_OSP(m_particles));
        }
        else if(boost::iequals(type, "PMB"))
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
        else if(boost::iequals(type, "gaussian bond force"))
        {
            forces.push_back(new PD_bondforceGaussian(m_particles));
        }
        else if(boost::iequals(type, "weighted bond"))
        {
            string weightType;
            cfg_forces[i].lookupValue("weightType", weightType);
            forces.push_back(new PD_bondforceGaussian(m_particles, weightType));
        }
        else if(boost::iequals(type, "DEM"))
        {
            double T, C;
            cfg_forces[i].lookupValue("T", T);
            cfg_forces[i].lookupValue("C", C);
            forces.push_back(new DemForce(m_particles, T, C));
        }
        else if(boost::iequals(type, "contact force"))
        {
            forces.push_back(new ContactForce(m_particles, m_grid, lc));
        }
        else if(boost::iequals(type, "viscous damper"))
        {
            double c;
            cfg_forces[i].lookupValue("c", c);
            forces.push_back(new ViscousDamper(m_particles, c));
        }
        else
        {
            cerr << "Force: " << type << " has not been implemented." << endl;
            exit(1);
        }

        forcesSet.push_back(type);
    }

    for(Force* force:forces)
    {
        vector<string> gParameters = force->initalGhostDependencies();
        for(const string &param:gParameters)
        {
            m_particles.addGhostParameter(param);
        }
    }

    if(isRoot)
        cout << "Initializing forces" << endl;

    // Initializing forces
    for(Force* force:forces)
    {
        force->numericalInitialization(calculateMicromodulus);
        force->initialize(E/E0, nu, delta, dim, h, lc);
    }


    if(isRoot)
    {
        cout << "Forces set: ";
        for(string f:forcesSet)
            cout << f << ", ";
        cout << endl;
    }
    //--------------------------------------------------------------------------
    // Recalcuating particle properties
    //--------------------------------------------------------------------------
    if(isRoot)
        cout << "Recalcuating particle properties" << endl;
    int applySurfaceCorrection = 0;
    m_cfg.lookupValue("applySurfaceCorrection", applySurfaceCorrection);
    if(applySurfaceCorrection)
    {
        for(Force* force:forces)
        {
            const int nSurfaceCorrections = 20;
            bool hasSurfaceCorrection = force->initializeSurfaceCorrection();

            if(!hasSurfaceCorrection)
                continue;
#if USE_MPI
            vector<string> additionGhostParameters = force->getSurfaceCorrectionGhostParameters();
            for(const string &param:additionGhostParameters)
            {
                m_particles.addGhostParameter(param);
            }

            m_grid.clearParticles();
            updateGrid(m_grid, m_particles);
            exchangeGhostParticles(m_grid, m_particles);
#endif
            for(int i=0; i<nSurfaceCorrections; i++)
            {
                force->applySurfaceCorrectionStep1(0.001);
#if USE_MPI
                m_grid.clearParticles();
                updateGrid(m_grid, m_particles);
                exchangeGhostParticles(m_grid, m_particles);
#endif
                force->applySurfaceCorrectionStep2();
            }
        }
    }

    if(!useS0fromCfg)
    {
        switch(dim)
        {
        case 3:
            reCalculatePdFractureCriterion(m_particles, G0, delta);
            break;
        case 2:
            reCalculatePdFractureCriterion(m_particles, G0, delta, h);
            break;
        default:
            cerr << "ERROR: dimension " << dim << " not supported for numerical 's0'" << endl;
            cerr << "use 2 or 3." << endl;
            exit(EXIT_FAILURE);
            break;
        }
    }

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

    if(boost::iequals(solverType, "ADR"))
    {
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
    else if(boost::iequals(solverType, "dynamic ADR"))
    {
        solver = new dynamicADR();
        m_cfg.lookupValue("dt", dt);
    }
    else if(boost::iequals(solverType, "conjugate gradient"))
    {
        const int maxIterations = 20000;
        const int threshold = 3.e-8;
        solver = new StaticSolver(maxIterations, threshold);
        m_cfg.lookupValue("dt", dt);
    }
    else if(boost::iequals(solverType, "velocity verlet"))
    {
        solver = new VelocityVerletIntegrator();
        m_cfg.lookupValue("dt", dt);
    }
    else if(boost::iequals(solverType, "euler-chromer"))
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
    solver->setDim(dim);
    solver->setRankAndCores(m_myRank, m_nCores);

    if(isRoot)
        cout << "Solver set: " << solverType << endl;
    //--------------------------------------------------------------------------
    // Setting the modifiers
    //--------------------------------------------------------------------------
    vector<Modifier *> modifiers;
    vector<Modifier *> spModifiers;
    vector<Modifier *> qsModifiers;
    vector<string> modifiersSet;


    Setting &cfg_modifiers = m_cfg.lookup("modifiers");
    for(int i=0; i<cfg_modifiers.getLength(); i++)
    {
        const char * tmpType;
        cfg_modifiers[i].lookupValue("type", tmpType);
        string type = tmpType;

        if(boost::iequals(type, "velocity Particles"))
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
        else if(boost::iequals(type, "move particles"))
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
        else if(boost::iequals(type, "boundary force"))
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

            if(boost::iequals(solverType, "ADR"))
            {
                modifiers.push_back(new boundaryForce(appliedForce, forceAxis, area, axis));
            }
            else
            {
                modifiers.push_back(new boundaryForce(appliedForce, forceAxis, area, axis));
            }
        }
        else if(boost::iequals(type, "PMB fracture"))
        {
            double alpha;

            if (!cfg_modifiers[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            if(boost::iequals(solverType, "ADR"))
            {
                qsModifiers.push_back(new ADRfracture(alpha));
            }
            else
            {
                spModifiers.push_back(new PmbFracture(alpha));
            }
        }
        else if(boost::iequals(type, "micropolar fracture"))
        {
            double criticalAngle;

            if (!cfg_modifiers[i].lookupValue("criticalAngle", criticalAngle))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            if(boost::iequals(solverType, "ADR"))
            {
                cerr << "Type: " << type << " not implemented for ADR" << endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                spModifiers.push_back(new MicropolarFracture(criticalAngle));
            }
        }
        else if(boost::iequals(type, "bond energy fracture"))
        {
            double G;

            if (!cfg_modifiers[i].lookupValue("G", G))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }

            Modifier * fMod = new BondEnergyFracture(delta, G, dim, forces, h);            
            if(boost::iequals(solverType, "ADR"))
            {
                qsModifiers.push_back(fMod);
            }
            else
            {
                spModifiers.push_back(fMod);
            }
        }
        else if(boost::iequals(type, "Mohr-Coulomb bond fracture"))
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
            if(boost::iequals(solverType, "ADR"))
            {
                AdrMohrCoulombBondFracture * failureCriterion = new AdrMohrCoulombBondFracture(mu, C, T, dim);
                qsModifiers.push_back(failureCriterion);
            }
            else
            {
                MohrCoulombBondFracture * failureCriterion = new MohrCoulombBondFracture(mu, C, T, dim);
                spModifiers.push_back(failureCriterion);
            }
        }
        else if(boost::iequals(type, "Mohr-Coulomb fracture"))
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
            if(boost::iequals(solverType, "ADR"))
            {
                ADRmohrCoulombFracture * failureCriterion = new ADRmohrCoulombFracture(mu, C, T, dim);
                qsModifiers.push_back(failureCriterion);
            }
            else
            {
                MohrCoulombFracture * failureCriterion = new MohrCoulombFracture(mu, C, T, dim);
                spModifiers.push_back(failureCriterion);
            }
        }
        else if(boost::iequals(type, "simple fracture"))
        {
            double alpha;

            if (!cfg_modifiers[i].lookupValue("alpha", alpha))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            SimpleFracture * simpleFracture = new SimpleFracture(alpha);            
            if(boost::iequals(solverType, "ADR"))
            {
                qsModifiers.push_back(simpleFracture);
            }
            else
            {
                spModifiers.push_back(simpleFracture);
            }
        }
        else if(boost::iequals(type, "ADR fracture average"))
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
        else if(boost::iequals(type, "rigid wall"))
        {
            int orientation;
            int topOrBottom;

            if (!cfg_modifiers[i].lookupValue("axis", orientation)
                || !cfg_modifiers[i].lookupValue("topOrBottom", topOrBottom))
            {
                cerr << "Error reading the parameters for modifier '"
                     << type << "'" << endl;
                exit(EXIT_FAILURE);
            }
            RigidWall * rigidWall = new RigidWall(m_grid, orientation, topOrBottom, lc);
            spModifiers.push_back(rigidWall);
        }
        else
        {
            cerr << "ERROR: modifier '" << type << "' does not exsits." << endl;
            exit(EXIT_FAILURE);
        }

        modifiersSet.push_back(type);
    }

    // Initializing modifiers
    for(Modifier* mod: modifiers)
    {
        mod->setDim(dim);
        mod->setParticles(m_particles);
        mod->registerParticleParameters();
        const auto & modNeededPorperties = mod->neededProperties();
        for(const auto & property:modNeededPorperties)
        {
            neededProperties.push_back(property);
        }

#if USE_MPI
        vector<string> additionGhostParameters = mod->initalGhostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
#endif
    }

    for(Modifier* mod: spModifiers)
    {
        mod->setDim(dim);
        mod->setParticles(m_particles);
        mod->registerParticleParameters();
        const auto & modNeededPorperties = mod->neededProperties();
        for(const auto & property:modNeededPorperties)
        {
            neededProperties.push_back(property);
        }
#if USE_MPI
        vector<string> additionGhostParameters = mod->initalGhostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
#endif
    }

    for(Modifier* mod: qsModifiers)
    {
        mod->setDim(dim);
        mod->setParticles(m_particles);
        mod->registerParticleParameters();
        const auto & modNeededPorperties = mod->neededProperties();
        for(const auto & property:modNeededPorperties)
        {
            neededProperties.push_back(property);
        }
#if USE_MPI
        vector<string> additionGhostParameters = mod->initalGhostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
#endif
    }
#if USE_MPI
    m_grid.clearGhostParticles();
    exchangeInitialGhostParticles(m_grid, m_particles);
#endif

    // Initializing modifiers
    for(Modifier* mod: modifiers)
    {
        mod->initialize();
    }

    for(Modifier* mod: spModifiers)
    {
        mod->initialize();
    }

    for(Modifier* mod: qsModifiers)
    {
        mod->initialize();
    }

    if(isRoot)
    {
        cout << "Modifiers set: ";
        for(string mod:modifiersSet)
            cout << mod << ", ";
        cout << endl;
    }

    //--------------------------------------------------------------------------
    // Setting additional initial conditions
    //--------------------------------------------------------------------------
    vector<string> intialConditionsSet;
    if( m_cfg.exists("initialConditions"))
    {
        Setting &cfg_initialConditions = m_cfg.lookup("initialConditions");
        for(int i=0; i<cfg_initialConditions.getLength(); i++)
        {
            const char * tmpType;
            cfg_initialConditions[i].lookupValue("type", tmpType);
            string type = tmpType;

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

            intialConditionsSet.push_back(type);
        }
    }

    if(isRoot)
    {
        cout << "Initial conditions set: ";
        for(string f:intialConditionsSet)
            cout << f << ", ";
        cout << endl;
    }
    //--------------------------------------------------------------------------
    // Setting the final ghost parameters
    //--------------------------------------------------------------------------
#if USE_MPI
    m_particles.clearGhostParameters();

    for(Force* force: forces)
    {
        vector<string> additionGhostParameters = force->ghostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
    }

    // Initializing modifiers
    for(Modifier* mod: modifiers)
    {
        vector<string> additionGhostParameters = mod->ghostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
    }

    for(Modifier* mod: spModifiers)
    {
        vector<string> additionGhostParameters = mod->ghostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
    }

    for(Modifier* mod: qsModifiers)
    {
        vector<string> additionGhostParameters = mod->ghostDependencies();
        for(const string &param:additionGhostParameters)
        {
            m_particles.addGhostParameter(param);
        }
    }

#endif
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
    int saveBinary = false;
    m_cfg.lookupValue("saveBinary", saveBinary);    

    SavePdData *saveParticles = new SavePdData(saveParameters);
    saveParticles->setDim(dim);
    saveParticles->setRankAndCores(m_myRank, m_nCores);
    saveParticles->setUpdateFrquency(saveFrequency);
    saveParticles->setGrid(&m_grid);
    saveParticles->setScaling(E0, L0, v0, t0, rho0);
    saveParticles->setSavePath(savePath);
    saveParticles->setParticles(&m_particles);
    saveParticles->setForces(forces);
    saveParticles->setWriteBinary(saveBinary);
    saveParticles->initialize();

    solver->setSaveInterval(saveFrequency);
    solver->setSaveParticles(saveParticles);

    const auto & saveNeededPorperties = saveParticles->neededProperties();
    for(const auto & property:saveNeededPorperties)
    {
        neededProperties.push_back(property);
    }
    //--------------------------------------------------------------------------
    // Setting the properties needs to be calculated
    //--------------------------------------------------------------------------
    vector<string> computeProperties;
    vector<CalculateProperty*> calcProperties;
    for(const auto prop:neededProperties)
    {
        const string type = prop.first;
        const int updateFrquency = prop.second;

        bool alreadyAdded = false;
        for(CalculateProperty* property:calcProperties)
        {
            if(boost::iequals(type, property->type))
            {
                if(property->updateFrquency() > updateFrquency)
                    property->setUpdateFrquency(updateFrquency);
                alreadyAdded = true;
                break;
            }
        }

        if(alreadyAdded)
            continue;

        if (boost::iequals(type, "PdAngle"))
        {
            CalculateProperty *property = new CalculatePdAngles();
            property->setUpdateFrquency(updateFrquency);
            calcProperties.push_back(property);
        }
        if(boost::iequals(type, "stress"))
        {
            CalculateProperty *property = new CalculateStress(forces);
            property->setUpdateFrquency(updateFrquency);
            calcProperties.push_back(property);
        }
        else
        {
            cerr << "Compute property '" << type << "' has not been implemented"
                 << endl;
            exit(1);
        }
        computeProperties.push_back(type);
    }

    for(auto prop:calcProperties)
    {
        prop->setDim(dim);
        prop->setParticles(m_particles);
        prop->initialize();
    }

    solver->setCalculateProperties(calcProperties);

    if(isRoot)
    {
        cout << "Compute properties: ";
        for(string f:computeProperties)
            cout << f << ", ";
        cout << endl;
    }


    //--------------------------------------------------------------------------
#if USE_MPI
    if(isRoot)
    {
        cout << "Ghost parameters: ";
        const vector<string> & gp =m_particles.ghostParametersString();
        for(const string &param:gp)
        {
            cout << param << ", ";
        }
        cout << endl;
    }
#endif
    //--------------------------------------------------------------------------
    double nSec = timer.toc();
    if(isRoot)
        cout << "Time: " << nSec << "s" << endl;



    return 0;
}
//------------------------------------------------------------------------------
void PdSolver::solve()
{
    if(isRoot)
        cout << "Starting solver" << endl;
    solver->solve();
}
//------------------------------------------------------------------------------
