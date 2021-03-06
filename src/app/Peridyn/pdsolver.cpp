#include "app/Peridyn/pdsolver.h"

#include "PDtools/SavePdData/savepddata.h"
#include <PDtools/CalculateProperties/calculateproperties.h>
#include <PDtools/Force/forces.h>
#include <PDtools/Modfiers/modifiers.h>
#include <PDtools/PdFunctions/pdfunctions.h>
#include <PDtools/Solver/solvers.h>

#include "Mesh/loadmesh.h"
#include "Mesh/meshtopdpartices.h"
#include "Mesh/pdmesh.h"

#include <boost/algorithm/string.hpp>
#include <boost/regex.h>

#ifdef USE_MPI
#include <mpi.h>
#endif


//------------------------------------------------------------------------------
PdSolver::PdSolver(std::string cfgPath, int myRank, int nMpiNodes)
    : m_configPath(cfgPath), m_myRank(myRank), m_nCores(nMpiNodes) {
  // Reading configuration file
  try {
    m_cfg.readFile(cfgPath.c_str());
  } catch (const libconfig::FileIOException &fioex) {
    std::cerr << "I/O error while reading the configuration file." << std::endl;
    exit(EXIT_FAILURE);
  } catch (const libconfig::ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    exit(EXIT_FAILURE);
  }
  m_cfg.setAutoConvert(true);
  if (myRank == 0) {
    isRoot = true;
  }
}
//------------------------------------------------------------------------------
PdSolver::~PdSolver() { delete solver; }
//------------------------------------------------------------------------------
int PdSolver::initialize() {
  // The entire program is initilized from the configuration file.
  // TODO(SBS): split the initialization into smaller functions.
  using namespace PDtools;
  arma::wall_clock timer;
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
  int planeStress = 1;
  int removeBondsOverVoids = 1;

  m_cfg.lookupValue("dim", dim);
  m_cfg.lookupValue("E", E);
  m_cfg.lookupValue("nu", nu);
  m_cfg.lookupValue("G0", G0);
  m_cfg.lookupValue("rho", rho);
  m_cfg.lookupValue("planeStress", planeStress);
  m_cfg.lookupValue("removeBondsOverVoids", removeBondsOverVoids);

  // Checking for dimensional scaling
  int scaleparameters = 0;
  m_cfg.lookupValue("scaleparameters", scaleparameters);

  if (scaleparameters) {
    m_cfg.lookupValue("E0", E0);
    m_cfg.lookupValue("L0", L0);
    m_cfg.lookupValue("v0", v0);
    m_cfg.lookupValue("t0", t0);
    m_cfg.lookupValue("rho0", rho0);
    E /= E0;
    rho /= rho0;
  }

  //--------------------------------------------------------------------------
  // Setting the domain
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Setting the domain" << endl;

  if (!m_cfg.exists("domain")) {
    cerr << "'domain' must be set in the configuration file" << endl;
    exit(EXIT_FAILURE);
  }
  libconfig::Setting &cfg_domain = m_cfg.lookup("domain");

  if (cfg_domain.getLength() != 6) {
    cerr << "'domain' must be set on the form ";
    cerr << "'domain = [x0, x1, y0, y1, z0, z1]'" << endl;
    exit(EXIT_FAILURE);
  }

  double dxdydz[M_DIM];
  vector<pair<double, double>> domain;
  for (int d = 0; d < M_DIM; d++) {
    pair<double, double> bound(cfg_domain[2 * d], cfg_domain[2 * d + 1]);
    bound.first /= L0;
    bound.second /= L0;
    dxdydz[d] = bound.second - bound.first;
    domain.push_back(bound);
  }

  // Setting periodicity
  arma::ivec3 periodicBoundaries = {0, 0, 0};

  try {
    libconfig::Setting &cfg_periodic = m_cfg.lookup("periodic");
    for (int d = 0; d < M_DIM; d++) {
      periodicBoundaries(d) = static_cast<int>(cfg_periodic[d]);
    }
  } catch (libconfig::SettingNotFoundException s) {
    if (m_myRank == 0)
      cout << "No periodicity set." << endl;
  }

  //--------------------------------------------------------------------------
  // Setting the grid
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Setting the grid" << endl;
  double delta;
  if (!m_cfg.lookupValue("delta", delta)) {
    cerr << "'delta' must be set in the configuration file" << endl;
    exit(EXIT_FAILURE);
  }
  double lc;
  if (!m_cfg.lookupValue("lc", lc)) {
    cerr << "'lc' must be set in the configuration file" << endl;
    exit(EXIT_FAILURE);
  }
  lc /= L0;
  delta /= L0;

  //    double gridspacing = 1.25*(delta + 0.5*lc);
  double gridspacing = 1.45 * (delta + 0.5 * lc);
  m_grid = Grid(domain, gridspacing, periodicBoundaries);
  m_grid.setIdAndCores(m_myRank, m_nCores);
  m_grid.dim(dim);
  m_grid.initialize();
  m_grid.setMyGridpoints();
  m_grid.setInitialPositionScaling(L0);
  if (m_myRank == 0) {
    vector<int> cpuConfig = m_grid.nCpuGrid();
    cout << "Running core configuration: ";
    for (int n : cpuConfig)
      cout << n << " ";
    cout << endl;
  }
  //--------------------------------------------------------------------------
  // Loading the particles
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Loading the particles" << endl;
  if (!m_cfg.exists("particlesPath")) {
    cerr << "'particlesPath' must be set in the configuration file" << endl;
    exit(EXIT_FAILURE);
  }
  string particlesPath =
      static_cast<const char *>(m_cfg.lookup("particlesPath"));
  //    string particlesPath = static_cast<const char
  //    *>(m_cfg.lookup("particlesPath"));
  string fileType = getFileEnding(particlesPath);

  if (boost::iequals(fileType, "msh")) {
    int quadratureDegree = 1;
    m_cfg.lookupValue("quadratureDegree", quadratureDegree);
    PdMesh msh = loadMesh2d(particlesPath);
    m_particles = convertMshToPdParticles(dim, quadratureDegree, msh, m_grid);
  } else {
    m_particles = load_pd(particlesPath, m_grid);
  }

  m_particles.dimensionalScaling(E0, L0, v0, t0, rho0);
  m_particles.dim(dim); // Brute forcing the dimension
  //--------------------------------------------------------------------------
  // TODO: Setting the initial position. Should not be done here
  //--------------------------------------------------------------------------
  mat &r0 = m_particles.r0();
  const mat &r = m_particles.r();

  for (unsigned int i = 0; i < m_particles.nParticles(); i++) {
    for (int d = 0; d < M_DIM; d++) {
      r0(i, d) = r(i, d);
    }
  }

  //--------------------------------------------------------------------------
  // Setting particles parameters
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Setting particles parameters" << endl;
  double h = dxdydz[2];
  double s0 = 1;
  bool useS0fromCfg = false;
  if (m_cfg.lookupValue("s0", s0)) {
    useS0fromCfg = true;
  }
  m_particles.registerParameter("rho", rho);
  m_particles.registerParameter("s0", s0);
  m_particles.registerParameter("radius");
  m_particles.registerParameter("groupId");

  calculateRadius(m_particles, dim, dxdydz[2]);

  int calculateMicromodulus = 0;
  m_cfg.lookupValue("calculateMicromodulus", calculateMicromodulus);

  //--------------------------------------------------------------------------
  // Setting the PD-connections
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Gridding particles and setting PD-connections" << endl;

  m_grid.clearParticles();
  m_grid.placeParticlesInGrid(m_particles);

// lc *= 1.05;
#ifdef USE_MPI
  vector<string> ghostParameters = {"volume", "radius", "groupId"};
  for (const string &param : ghostParameters) {
    m_particles.addGhostParameter(param);
  }
  exchangeGhostParticles(m_grid, m_particles);
#endif
  int performVolumeCorrection = 1;
  m_cfg.lookupValue("performVolumeCorrection", performVolumeCorrection);
  setPdConnections(m_particles, m_grid, delta, lc);
  m_grid.clearGhostParticles();
#if USE_MPI
  exchangeInitialGhostParticles(m_grid, m_particles);
#endif

  if (removeBondsOverVoids) {
    m_grid.clearParticles();
    updateGrid(m_grid, m_particles, true);
#if USE_MPI
    exchangeInitialGhostParticles(m_grid, m_particles);
#endif
    removeVoidConnections(m_particles, m_grid, delta, lc);
  }
  cleanUpPdConnections(m_particles);

  m_particles.registerPdParameter("volumeScaling", 1);
  if (performVolumeCorrection) {
    applyVolumeCorrection(m_particles, delta, lc, dim);
  }

  setPD_N3L(m_particles);
  //--------------------------------------------------------------------------
  // Setting the solver
  //--------------------------------------------------------------------------
  string solverType = static_cast<const char *>(m_cfg.lookup("solverType"));

  int nSteps;
  double dt;

  if (!m_cfg.lookupValue("nSteps", nSteps)) {
    cerr << "Error reading the 'nSteps' in config file" << endl;
    exit(EXIT_FAILURE);
  }

  if (boost::iequals(solverType, "ADR")) {
    ADR *adrSolver = new ADR();
    dt = 1.0;
    double errorThreshold;
    int maxSteps = 5000;
    int maxStepsFracture = 1000;
    m_cfg.lookupValue("maxSteps", maxSteps);
    m_cfg.lookupValue("maxStepsFracture", maxStepsFracture);

    if (!m_cfg.lookupValue("errorThreshold", errorThreshold)) {
      cerr << "Error reading the 'errorThreshold' in config file" << endl;
      exit(EXIT_FAILURE);
    }

    adrSolver->maxSteps(maxSteps);
    adrSolver->maxStepsFracture(maxStepsFracture);
    adrSolver->setErrorThreshold(errorThreshold);
    solver = adrSolver;
  } else if (boost::iequals(solverType, "dynamic ADR")) {
    solver = new dynamicADR();
    m_cfg.lookupValue("dt", dt);
  } else if (boost::iequals(solverType, "conjugate gradient")) {
    const int maxIterations = 20000;
    const int threshold = 3.e-8;
    solver = new StaticSolver(maxIterations, threshold);
    m_cfg.lookupValue("dt", dt);
  } else if (boost::iequals(solverType, "velocity verlet")) {
    solver = new VelocityVerletIntegrator();
    m_cfg.lookupValue("dt", dt);
  } else if (boost::iequals(solverType, "euler-chromer")) {
    solver = new EulerCromerIntegrator();
    m_cfg.lookupValue("dt", dt);
  } else {
    cerr << "Error: solver not set" << endl;
    exit(EXIT_FAILURE);
  }
  dt /= t0;
  solver->setDim(dim);
  solver->setRankAndCores(m_myRank, m_nCores);

  if (isRoot)
    cout << "Solver set: " << solverType << endl;

  //--------------------------------------------------------------------------
  // Setting the Forces
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Setting the Forces" << endl;

  vector<string> forcesSet;
  vector<Force *> forces;
  libconfig::Setting &cfg_forces = m_cfg.lookup("forces");

  for (int i = 0; i < cfg_forces.getLength(); i++) {
    const char *tmpType;
    cfg_forces[i].lookupValue("type", tmpType);
    string type = tmpType;

    // -------- BOND MODELS --------
    if (boost::iequals(type, "bond force")) {
      forces.push_back(new PD_bondForce(m_particles));
    } else if (boost::iequals(type, "dampened bond force")) {
      double c;
      if (!cfg_forces[i].lookupValue("c", c)) {
        cerr << "The dampening parameter 'c' must be set for " << type << endl;
        exit(EXIT_FAILURE);
      }
      forces.push_back(new PD_dampenedBondForce(m_particles, c / t0));
    } else if (boost::iequals(type, "PMB")) {
      double alpha;

      if (!cfg_forces[i].lookupValue("alpha", alpha)) {
        cerr << "Error reading the parameters for force '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      forces.push_back(new PD_PMB(m_particles, lc, delta, alpha));
    } else if (boost::iequals(type, "PMB LIN")) {
      double alpha;

      if (!cfg_forces[i].lookupValue("alpha", alpha)) {
        cerr << "Error reading the parameters for force '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      forces.push_back(
          new PD_PMB_LINEAR_INTEGRATOR(m_particles, lc, delta, alpha));
    } else if (boost::iequals(type, "gaussian bond force")) {
      forces.push_back(new PD_bondforceGaussian(m_particles));
    } else if (boost::iequals(type, "weighted bond")) {
      string weightType;
      cfg_forces[i].lookupValue("weightType", weightType);
      forces.push_back(new PD_bondforceGaussian(m_particles, weightType));
    }
    // -------- LPS MODELS --------
    else if (boost::iequals(type, "LPS")) {
      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);
      //      forces.push_back(new PD_LPS2(m_particles, planeStress,
      //      analyticalM));
      forces.push_back(new PD_LPS(m_particles, planeStress, analyticalM));
    } else if (boost::iequals(type, "LPS SHEAR")) {
      double G0;
      if (!cfg_forces[i].lookupValue("G0", G0)) {
        int useOptimized = false;

        if (!cfg_forces[i].lookupValue("useOptimized", useOptimized))
          forces.push_back(new PD_LPSS(m_particles, planeStress));
        else
          forces.push_back(new PD_LPSS_opt(m_particles, planeStress));
      } else {
        forces.push_back(new PD_LPSS_G(m_particles, G0, planeStress));
      }
    } else if (boost::iequals(type, "LPS K")) {
      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);
      forces.push_back(new PD_LPS_K(m_particles, planeStress, analyticalM));
    } else if (boost::iequals(type, "dampened LPS")) {
      double c;
      if (!cfg_forces[i].lookupValue("c", c)) {
        cerr << "The dampening parameter 'c' must be set for " << type << endl;
        exit(EXIT_FAILURE);
      }
      forces.push_back(
          new PD_lpsDampenedContact(m_particles, c / t0, planeStress));
    } else if (boost::iequals(type, "LPS model")) {
      double c, mu, C, T;

      if ((!cfg_forces[i].lookupValue("c", c) &&
           !boost::iequals(solverType, "ADR")) ||
          !cfg_forces[i].lookupValue("mu", mu) ||
          !cfg_forces[i].lookupValue("C", C) ||
          !cfg_forces[i].lookupValue("T", T)) {
        cerr << "The dampening parameter 'c' must be set for " << type << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;

      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);

      if (boost::iequals(solverType, "ADR")) {
        forces.push_back(
            new PD_LPS_adrmc(m_particles, mu, C, T, planeStress, analyticalM));
      } else {
        forces.push_back(new LPS_mc(m_particles, c / t0, mu, C, T, planeStress,
                                    analyticalM));
      }
    } else if (boost::iequals(type, "LPS stretch")) {
      double c, stretchCrit, shearCrit;

      if ((!cfg_forces[i].lookupValue("c", c) &&
           !boost::iequals(solverType, "ADR")) ||
          !cfg_forces[i].lookupValue("stretchCrit", stretchCrit) ||
          !cfg_forces[i].lookupValue("shearCrit", shearCrit)) {
        cerr << "The dampening parameter 'c' must be set for " << type << endl;
        exit(EXIT_FAILURE);
      }
      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);

      if (boost::iequals(solverType, "ADR")) {
        forces.push_back(new PD_LPS_ADR_STRAIN(m_particles, c / t0, stretchCrit,
                                               shearCrit, planeStress,
                                               analyticalM));
      } else {
        forces.push_back(new PD_LPS_CRIT_STRAIN(m_particles, c / t0,
                                                stretchCrit, shearCrit,
                                                planeStress, analyticalM));
      }
    } // -------- NOSBPD MODELS --------
    else if (boost::iequals(type, "NOSBPD MC")) {
      double mu, C, T;
      if (!cfg_forces[i].lookupValue("mu", mu) ||
          !cfg_forces[i].lookupValue("C", C) ||
          !cfg_forces[i].lookupValue("T", T)) {
        cerr << "The dampening parameter 'c' must be set for " << type << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;

      forces.push_back(new PD_NOPD(m_particles, mu, C, T, planeStress));
    } else if (boost::iequals(type, "OSP")) {
      forces.push_back(new PD_OSP(m_particles));
    }
    // -------- OTHER FORCES --------
    else if (boost::iequals(type, "DEM")) {
      double T, C;
      cfg_forces[i].lookupValue("T", T);
      cfg_forces[i].lookupValue("C", C);
      forces.push_back(new DemForce(m_particles, T, C));
    } else if (boost::iequals(type, "contact force")) {
      //      double interactionRadius = 0.95;
      //      double interactionScaling = 15;
      int updateFrquency = 10;
      cfg_forces[i].lookupValue("verletUpdateFrq", updateFrquency);
      forces.push_back(
          new ContactForce(m_particles, m_grid, lc, updateFrquency));
    } else if (boost::iequals(type, "viscous damper")) {
      double c;
      cfg_forces[i].lookupValue("c", c);
      forces.push_back(new ViscousDamper(m_particles, c));
    }
    // -------- LPS POROSITY MODELS --------
    else if (boost::iequals(type, "LPS porosity")) {
      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);
      double phi_c, n;
      if (!cfg_forces[i].lookupValue("phi_c", phi_c) ||
          !cfg_forces[i].lookupValue("n", n)) {
        cerr << "The material parameter 'phi_c' and 'n' must be set"
             << "They are used to scale material parameters. Example:"
             << "E = E0(1-phi/phi_c)^n" << endl
             << "nu = nu0(1-phi/phi_c)." << endl
             << type << endl;
        exit(EXIT_FAILURE);
      }

      forces.push_back(
          new PD_LPS_POROSITY(m_particles, phi_c, n, planeStress, analyticalM));
    } else if (boost::iequals(type, "dampened LPS porosity")) {
      double c;
      double phi_c, n;

      if (!cfg_forces[i].lookupValue("c", c) ||
          !cfg_forces[i].lookupValue("phi_c", phi_c) ||
          !cfg_forces[i].lookupValue("n", n)) {
        cerr << "The dampening parameter 'c' must be set."
             << "The material parameter 'phi_c' and 'n' must be set"
             << "They are used to scale material parameters. Example:"
             << "E = E0(1-phi/phi_c)^n" << endl
             << "nu = nu0(1-phi/phi_c)." << endl
             << type << endl;
        exit(EXIT_FAILURE);
      }
      forces.push_back(new PD_lpsDampenedContact_porosity(m_particles, phi_c, n,
                                                          c / t0, planeStress));
    } else if (boost::iequals(type, "LPS model porosity")) {
      double c, mu, C, T;
      double m, b;

      if ((!cfg_forces[i].lookupValue("c", c) &&
           !boost::iequals(solverType, "ADR")) ||
          !cfg_forces[i].lookupValue("mu", mu) ||
          !cfg_forces[i].lookupValue("C", C) ||
          !cfg_forces[i].lookupValue("T", T) ||
          !cfg_forces[i].lookupValue("m", m) ||
          !cfg_forces[i].lookupValue("b", b)) {
        cerr << "The dampening parameter 'c' must be set " << endl;
        cerr << "The fracture parameters 'mu', 'C' and 'T' must be set for "
             << endl;
        cerr << "The porosity parameters 'm' and 'b' must be set.";
        cerr << "They are used to scale material parameters. Example:"
             << " E = E0*exp(b(m-phi))." << endl
             << type << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;

      int analyticalM = false;
      cfg_forces[i].lookupValue("analyticalM", analyticalM);

      if (boost::iequals(solverType, "ADR")) {
        forces.push_back(new PD_LPS_porosity_adrmc(m_particles, m, b, mu, C, T,
                                                   planeStress, analyticalM));
      } else {
        forces.push_back(new LPS_porosity_mc(m_particles, m, b, c / t0, mu, C,
                                             T, planeStress, analyticalM));
      }
    }

    else {
      cerr << "Force: " << type << " has not been implemented." << endl;
      exit(1);
    }

    forcesSet.push_back(type);
  }

  for (Force *force : forces) {
    vector<string> gParameters = force->initalGhostDependencies();
    for (const string &param : gParameters) {
      m_particles.addGhostParameter(param);
    }
  }

  if (!useS0fromCfg) {
    /*
    switch(dim)
    {
    case 3:
        reCalculatePdFractureCriterion(m_particles, G0, delta);
        break;
    case 2:
        reCalculatePdFractureCriterion(m_particles, G0, delta, h);
        break;
    default:
        cerr << "ERROR: dimension " << dim << " not supported for numerical
    's0'" << endl;
        cerr << "use 2 or 3." << endl;
        exit(EXIT_FAILURE);
        break;
    }
    */
  }

  //--------------------------------------------------------------------------
  // Setting the modifiers
  //--------------------------------------------------------------------------
  //    MPI_Barrier(MPI_COMM_WORLD);
  if (isRoot) {
    cout << "Setting modifiers:" << endl;
  }
  vector<Modifier *> boundaryModifiers;
  vector<Modifier *> spModifiers;
  vector<Modifier *> qsModifiers;
  vector<string> modifiersSet;

  libconfig::Setting &cfg_modifiers = m_cfg.lookup("modifiers");
  for (int i = 0; i < cfg_modifiers.getLength(); i++) {
    const char *tmpType;
    cfg_modifiers[i].lookupValue("type", tmpType);
    string type = tmpType;

    if (boost::iequals(type, "velocity Particles")) {
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
          !cfg_modifiers[i].lookupValue("vAxis", vAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      v /= v0;
      boundaryModifiers.push_back(
          new VelocityBoundary(v, vAxis, area, axis, dt, nSteps, isStatic));
    } else if (boost::iequals(type, "move particles")) {
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
          !cfg_modifiers[i].lookupValue("vAxis", vAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      v /= v0;
      boundaryModifiers.push_back(
          new MoveParticles(v, vAxis, area, axis, dt, isStatic));
    } else if (boost::iequals(type, "move particles zone")) {
      double v;
      int isStatic = 0;
      vector<double> area = {0, 0, 0, 0, 0, 0};
      vec velocityDirection = {0, 0, 0};
      for (int d = 0; d < dim; d++) {
        area[2 * d] = cfg_modifiers[i]["area"][2 * d];
        area[2 * d + 1] = cfg_modifiers[i]["area"][2 * d + 1];
        area[2 * d] /= L0;
        area[2 * d + 1] /= L0;
        velocityDirection(d) = cfg_modifiers[i]["velocityDirection"][d];
      }
      cfg_modifiers[i].lookupValue("static", isStatic);

      if (!cfg_modifiers[i].lookupValue("v", v)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      v /= v0;

      boundaryModifiers.push_back(new MoveParticlesZone(
          v, velocityDirection, area, dt, isStatic, delta));
    } else if (boost::iequals(type, "move particle group")) {
      int groupId;
      double v;
      int isStatic = 0;
      vec velocityDirection = {0, 0, 0};
      for (int d = 0; d < dim; d++) {
        velocityDirection(d) = cfg_modifiers[i]["velocityDirection"][d];
      }
      cfg_modifiers[i].lookupValue("static", isStatic);

      if (!cfg_modifiers[i].lookupValue("v", v) ||
          !cfg_modifiers[i].lookupValue("groupId", groupId)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      v /= v0;

      boundaryModifiers.push_back(
          new MoveParticleGroup(groupId, v, velocityDirection, dt, isStatic));
    } else if (boost::iequals(type, "strain particles")) {
      double strainrate;
      int axis, strainAxis;
      int isStatic = 0;

      double a0 = cfg_modifiers[i]["area"][0];
      double a1 = cfg_modifiers[i]["area"][1];
      cfg_modifiers[i].lookupValue("static", isStatic);
      a0 /= L0;
      a1 /= L0;

      pair<double, double> area(a0, a1);

      if (!cfg_modifiers[i].lookupValue("strainrate", strainrate) ||
          !cfg_modifiers[i].lookupValue("axis", axis) ||
          !cfg_modifiers[i].lookupValue("strainAxis", strainAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }

      boundaryModifiers.push_back(new StrainBoundary(strainrate, strainAxis, nu,
                                                     area, axis, dt, isStatic));
    } else if (boost::iequals(type, "boundary force")) {
      double appliedForce;
      int axis, forceAxis;
      int incremental = 1;
      int steps = -1;
      double a0 = cfg_modifiers[i]["area"][0];
      double a1 = cfg_modifiers[i]["area"][1];
      a0 /= L0;
      a1 /= L0;

      pair<double, double> area(a0, a1);
      cfg_modifiers[i].lookupValue("steps", steps);
      cfg_modifiers[i].lookupValue("incremental", incremental);

      if (!cfg_modifiers[i].lookupValue("appliedForce", appliedForce) ||
          !cfg_modifiers[i].lookupValue("axis", axis) ||
          !cfg_modifiers[i].lookupValue("forceAxis", forceAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      appliedForce /= (E0 / pow(L0, 4));
      cout << "incremental:" << incremental << endl;

      double volume =
          (a1 - a0) * dxdydz[0] * dxdydz[1] * dxdydz[2] / dxdydz[axis];
      appliedForce /= volume;
      boundaryModifiers.push_back(new boundaryForce(
          appliedForce, forceAxis, area, axis, steps, delta, incremental));
    } else if (boost::iequals(type, "boundary force scaled")) {
      double appliedForce;
      int axis, forceAxis;
      int incremental = 1;
      int scaled = 1.;
      int steps = -1;
      double a0 = cfg_modifiers[i]["area"][0];
      double a1 = cfg_modifiers[i]["area"][1];
      a0 /= L0;
      a1 /= L0;

      pair<double, double> area(a0, a1);
      cfg_modifiers[i].lookupValue("steps", steps);
      cfg_modifiers[i].lookupValue("incremental", incremental);

      if (!cfg_modifiers[i].lookupValue("appliedForce", appliedForce) ||
          !cfg_modifiers[i].lookupValue("axis", axis) ||
          !cfg_modifiers[i].lookupValue("forceAxis", forceAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      appliedForce /= (E0 / pow(L0, 4));
      cout << "incremental:" << incremental << endl;

      boundaryModifiers.push_back(new boundaryForce(appliedForce, forceAxis,
                                                    area, axis, steps, delta,
                                                    incremental, scaled));

    } else if (boost::iequals(type, "boundary stress")) {
      double appliedStress;
      int axis, stressAxis;
      int incremental = 1;
      int steps = -1;
      double a0 = cfg_modifiers[i]["area"][0];
      double a1 = cfg_modifiers[i]["area"][1];
      a0 /= L0;
      a1 /= L0;

      pair<double, double> area(a0, a1);
      cfg_modifiers[i].lookupValue("steps", steps);
      cfg_modifiers[i].lookupValue("incremental", incremental);

      if (!cfg_modifiers[i].lookupValue("stress", appliedStress) ||
          !cfg_modifiers[i].lookupValue("axis", axis) ||
          !cfg_modifiers[i].lookupValue("stressAxis", stressAxis)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }

      //            appliedStress /= (E0/pow(L0, 4));
      //            appliedStress /= (a1 - a0);
      appliedStress /= delta;
      boundaryModifiers.push_back(new boundaryForce(
          appliedStress, stressAxis, area, axis, steps, delta, incremental));
    } else if (boost::iequals(type, "PMB fracture")) {
      double alpha;

      if (!cfg_modifiers[i].lookupValue("alpha", alpha)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }

      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(new ADRfracture(alpha));
      } else {
        spModifiers.push_back(new PmbFracture(alpha));
      }
    } else if (boost::iequals(type, "bond energy fracture")) {
      double G;

      if (!cfg_modifiers[i].lookupValue("G", G)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }

      Modifier *fMod = new BondEnergyFracture(delta, G, forces, h);
      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(fMod);
      } else {
        spModifiers.push_back(fMod);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb bond fracture")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      if (boost::iequals(solverType, "ADR")) {
        AdrMohrCoulombBondFracture *failureCriterion =
            new AdrMohrCoulombBondFracture(mu, C, T);
        qsModifiers.push_back(failureCriterion);
      } else {
        MohrCoulombBondFracture *failureCriterion =
            new MohrCoulombBondFracture(mu, C, T);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb fracture")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      if (boost::iequals(solverType, "ADR")) {
        ADRmohrCoulombFracture *failureCriterion =
            new ADRmohrCoulombFracture(mu, C, T);
        qsModifiers.push_back(failureCriterion);
      } else {
        MohrCoulombFracture *failureCriterion =
            new MohrCoulombFracture(mu, C, T);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb node split")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      if (boost::iequals(solverType, "ADR")) {
        cerr << "ADR not implemented '" << type << "'" << endl;
        exit(EXIT_FAILURE);
      } else {
        MohrCoulombNodeSplit *failureCriterion =
            new MohrCoulombNodeSplit(mu, C, T);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb max connected")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      if (boost::iequals(solverType, "ADR")) {
        cerr << "ADR not implemented '" << type << "'" << endl;
        exit(EXIT_FAILURE);
      } else {
        MohrCoulomMaxConnected *failureCriterion =
            new MohrCoulomMaxConnected(mu, C, T);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb weighted average")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;

      if (boost::iequals(solverType, "ADR")) {
        //                ADRmohrCoulombFracture * failureCriterion = new
        //                ADRmohrCoulombFracture(mu, C, T, dim);
        //                qsModifiers.push_back(failureCriterion);
      } else {
        MohrCoulombWeightedAverage *failureCriterion =
            new MohrCoulombWeightedAverage(mu, C, T);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb max")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      MohrCoulombMax *failureCriterion = new MohrCoulombMax(mu, C, T);

      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(failureCriterion);
      } else {
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb max fracture")) {
      double mu, C, T;

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;
      MohrCoulombMaxFracture *failureCriterion =
          new MohrCoulombMaxFracture(mu, C, T);

      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(failureCriterion);
      } else {
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Mohr-Coulomb max fracture weighted")) {
      double mu, C, T;
      double Wc = 0.9;
      double Bc = 0.9;
      cfg_modifiers[i].lookupValue("W", Wc);
      cfg_modifiers[i].lookupValue("B", Bc);

      if (!cfg_modifiers[i].lookupValue("mu", mu) ||
          !cfg_modifiers[i].lookupValue("C", C) ||
          !cfg_modifiers[i].lookupValue("T", T)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      C /= E0;
      T /= E0;

      if (boost::iequals(solverType, "ADR")) {
        MohrCoulombMaxFractureWeighted *failureCriterion =
            new MohrCoulombMaxFractureWeightedAdr(mu, C, T, Wc, Bc);
        qsModifiers.push_back(failureCriterion);
      } else {
        MohrCoulombMaxFractureWeighted *failureCriterion =
            new MohrCoulombMaxFractureWeighted(mu, C, T, Wc, Bc);
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "Strain fracture")) {
      double Eeq, Evol;

      if (!cfg_modifiers[i].lookupValue("Eeq", Eeq) ||
          !cfg_modifiers[i].lookupValue("Evol", Evol)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }

      StrainFracture *failureCriterion = new StrainFracture(Eeq, Evol, dim);

      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(failureCriterion);
      } else {
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "simple fracture")) {
      double alpha;

      if (!cfg_modifiers[i].lookupValue("alpha", alpha)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      SimpleFracture *failureCriterion = new SimpleFracture(alpha);
      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(failureCriterion);
      } else {
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "von Mises fracture")) {
      double sigma_y;

      if (!cfg_modifiers[i].lookupValue("sigma_y", sigma_y)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      VonMisesFracture *failureCriterion = new VonMisesFracture(sigma_y);
      if (boost::iequals(solverType, "ADR")) {
        qsModifiers.push_back(failureCriterion);
      } else {
        spModifiers.push_back(failureCriterion);
      }
    } else if (boost::iequals(type, "ADR fracture average")) {
      double alpha;

      if (!cfg_modifiers[i].lookupValue("alpha", alpha)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      ADRfractureAverage *adrFracture = new ADRfractureAverage(alpha);
      qsModifiers.push_back(adrFracture);
    } else if (boost::iequals(type, "rigid wall")) {
      int orientation;
      int topOrBottom;

      if (!cfg_modifiers[i].lookupValue("axis", orientation) ||
          !cfg_modifiers[i].lookupValue("topOrBottom", topOrBottom)) {
        cerr << "Error reading the parameters for modifier '" << type << "'"
             << endl;
        exit(EXIT_FAILURE);
      }
      RigidWall *rigidWall =
          new RigidWall(m_grid, orientation, topOrBottom, lc);
      spModifiers.push_back(rigidWall);
    } else if (boost::iequals(type, "custom fracture")) {
      //            vector<double> fracture = {0., 0.3, 0.09999, 0.100001}; //
      //            Test fracture
      vector<double> fracture = {0., 0., 0.,
                                 0.}; // Test fracture, [0, x1, y0, y1], all z

      fracture[0] = cfg_modifiers[i]["fracture"][0];
      fracture[1] = cfg_modifiers[i]["fracture"][1];
      fracture[2] = cfg_modifiers[i]["fracture"][2];
      fracture[3] = cfg_modifiers[i]["fracture"][3];
      addFractures(m_particles, domain, fracture);
    } else {
      cerr << "ERROR: modifier '" << type << "' does not exsits." << endl;
      exit(EXIT_FAILURE);
    }

    modifiersSet.push_back(type);
  }
  //    MPI_Barrier(MPI_COMM_WORLD);
  if (m_myRank == 0)
    cout << "Initializing modifiers " << endl;

  // Initializing modifiers
  for (Modifier *mod : boundaryModifiers) {
    mod->setDim(dim);
    mod->setGrid(&m_grid);
    mod->setParticles(m_particles);
    mod->registerParticleParameters();

    const auto &modNeededPorperties = mod->neededProperties();
    for (const auto &property : modNeededPorperties) {
      neededProperties.push_back(property);
    }
#if USE_MPI
    vector<string> additionGhostParameters = mod->initalGhostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
#endif
  }

  for (Modifier *mod : spModifiers) {
    mod->setDim(dim);
    mod->setGrid(&m_grid);
    mod->setParticles(m_particles);
    mod->registerParticleParameters();
    const auto &modNeededPorperties = mod->neededProperties();
    for (const auto &property : modNeededPorperties) {
      neededProperties.push_back(property);
    }
#if USE_MPI
    vector<string> additionGhostParameters = mod->initalGhostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
#endif
  }

  for (Modifier *mod : qsModifiers) {
    mod->setDim(dim);
    mod->setGrid(&m_grid);
    mod->setParticles(m_particles);
    mod->registerParticleParameters();
    const auto &modNeededPorperties = mod->neededProperties();
    for (const auto &property : modNeededPorperties) {
      neededProperties.push_back(property);
    }
#if USE_MPI
    vector<string> additionGhostParameters = mod->initalGhostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
#endif
  }
#if USE_MPI
  m_grid.clearGhostParticles();
  MPI_Barrier(MPI_COMM_WORLD);
  if (m_myRank == 0)
    cout << "exchangeInitialGhostParticles" << endl;
  exchangeInitialGhostParticles(m_grid, m_particles);
#endif

  //    MPI_Barrier(MPI_COMM_WORLD);
  //    if(m_myRank == 0)
  //       cout << " modifiers mod->initialize()" << endl;
  // Initializing modifiers
  for (Modifier *mod : boundaryModifiers) {
    mod->initialize();
  }

  for (Modifier *mod : spModifiers) {
    mod->initialize();
  }

  for (Modifier *mod : qsModifiers) {
    mod->initialize();
  }

  if (isRoot) {
    cout << "Modifiers set: ";
    for (string mod : modifiersSet)
      cout << mod << ", ";
    cout << endl;
  }

  //--------------------------------------------------------------------------
  // Initializing forces - moved here for fracture modifiers
  //--------------------------------------------------------------------------
  if (isRoot)
    cout << "Initializing forces" << endl;

  // Initializing forces
  for (Force *force : forces) {
    force->numericalInitialization(calculateMicromodulus);
    force->initialize(E, nu, delta, dim, h, lc);
    const auto &needed = force->getNeededProperties();

    for (auto prop : needed) {
      neededProperties.push_back(prop);
    }
  }

  if (isRoot) {
    cout << "Forces set: ";
    for (string f : forcesSet)
      cout << f << ", ";
    cout << endl;
  }
  //--------------------------------------------------------------------------
  // Setting additional initial conditions
  //--------------------------------------------------------------------------
  vector<string> intialConditionsSet;
  if (m_cfg.exists("initialConditions")) {
    libconfig::Setting &cfg_initialConditions =
        m_cfg.lookup("initialConditions");

    for (int i = 0; i < cfg_initialConditions.getLength(); i++) {
      const char *tmpType;
      cfg_initialConditions[i].lookupValue("type", tmpType);
      string type = tmpType;

      if (type == "strain") {
        double strain;
        int axis;

        double a0 = cfg_initialConditions[i]["area"][0];
        double a1 = cfg_initialConditions[i]["area"][1];
        a0 /= L0;
        a1 /= L0;

        pair<double, double> area(a0, a1);

        if (!cfg_initialConditions[i].lookupValue("strain", strain) ||
            !cfg_initialConditions[i].lookupValue("axis", axis) ||
            !cfg_initialConditions[i].lookupValue("axis", axis)) {
          cerr << "Error reading the parameters for modifier '" << type << "'"
               << endl;
          exit(EXIT_FAILURE);
        }
        applyInitialStrainStrain(m_particles, strain, axis, area);
      }

      intialConditionsSet.push_back(type);
    }
  }

  if (isRoot) {
    cout << "Initial conditions set: ";
    for (string f : intialConditionsSet)
      cout << f << ", ";
    cout << endl;
  }

//--------------------------------------------------------------------------
// Setting the final ghost parameters
//--------------------------------------------------------------------------
#if USE_MPI
  m_particles.clearGhostParameters();

  for (Force *force : forces) {
    vector<string> additionGhostParameters = force->ghostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
  }

  // Initializing modifiers
  for (Modifier *mod : boundaryModifiers) {
    vector<string> additionGhostParameters = mod->ghostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
  }

  for (Modifier *mod : spModifiers) {
    vector<string> additionGhostParameters = mod->ghostDependencies();
    for (const string &param : additionGhostParameters) {
      m_particles.addGhostParameter(param);
    }
  }

  for (Modifier *mod : qsModifiers) {
    vector<string> additionGhostParameters = mod->ghostDependencies();
    for (const string &param : additionGhostParameters) {
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

  for (Force *force : forces) {
    solver->addForce(force);
  }
  for (Modifier *mod : boundaryModifiers) {
    solver->addBoundaryModifier(mod);
  }
  for (Modifier *mod : spModifiers) {
    solver->addSpModifier(mod);
  }
  for (Modifier *mod : qsModifiers) {
    solver->addQsModifiers(mod);
  }

  //--------------------------------------------------------------------------
  // Setting the save values
  //--------------------------------------------------------------------------
  int saveFrequency;
  if (!m_cfg.lookupValue("saveFrequency", saveFrequency)) {
    cerr << "Error reading the 'saveFrequency' in config file" << endl;
    exit(EXIT_FAILURE);
  }
  if (!m_cfg.exists("savePath")) {
    cerr << "Error reading the 'savePath' in config file" << endl;
    exit(EXIT_FAILURE);
  }
  string savePath = static_cast<const char *>(m_cfg.lookup("savePath"));

  vector<string> saveParameters;
  if (m_cfg.exists("saveParameters")) {
    const libconfig::Setting &cfg_saveParameters =
        m_cfg.lookup("saveParameters");
    for (int i = 0; i < cfg_saveParameters.getLength(); i++) {
      saveParameters.push_back(cfg_saveParameters[i].c_str());
    }
  } else {
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

  const auto &saveNeededProperties = saveParticles->neededProperties();
  for (const auto &property : saveNeededProperties) {
    neededProperties.push_back(property);
  }

  //--------------------------------------------------------------------------
  // Setting the properties needs to be calculated
  //--------------------------------------------------------------------------
  vector<string> computeProperties;
  vector<CalculateProperty *> calcProperties;

  for (const auto prop : neededProperties) {
    string type = prop.first;
    const int updateFrquency = prop.second;

    bool alreadyAdded = false;

    if (boost::iequals(type, "strain")) {
      type = "stress";
    }

    for (CalculateProperty *property : calcProperties) {
      if (boost::iequals(type, property->type)) {
        if (property->updateFrequency() > updateFrquency)
          property->setUpdateFrquency(updateFrquency);
        alreadyAdded = true;
        break;
      }
    }

    if (alreadyAdded)
      continue;

    if (boost::iequals(type, "PdAngle")) {
      CalculateProperty *property = new CalculatePdAngles();
      property->setUpdateFrquency(updateFrquency);
      calcProperties.push_back(property);
    } else if (boost::iequals(type, "stress") ||
               boost::iequals(type, "strain")) {
      CalculateProperty *property =
          new CalculateStressStrain(forces, E, nu, delta, planeStress);
      property->setUpdateFrquency(updateFrquency);
      calcProperties.push_back(property);
    } else if (boost::iequals(type, "stress2")) {
      CalculateProperty *property = new CalculateStress(forces);
      property->setUpdateFrquency(updateFrquency);
      calcProperties.push_back(property);
    } else if (boost::iequals(type, "damage")) {
      CalculateProperty *property = new CalculateDamage(delta);
      property->setUpdateFrquency(updateFrquency);
      calcProperties.push_back(property);
    } else {
      cerr << "Compute property '" << type << "' has not been implemented"
           << endl;
      exit(1);
    }
    computeProperties.push_back(type);
  }

  for (auto prop : calcProperties) {
    prop->setDim(dim);
    prop->setParticles(m_particles);
    prop->initialize();
  }
  solver->setCalculateProperties(calcProperties);
  if (isRoot) {
    cout << "Compute properties: ";
    for (string f : computeProperties)
      cout << f << ", ";
    cout << endl;
  }
//--------------------------------------------------------------------------
#if USE_MPI
  if (isRoot) {
    cout << "Ghost parameters: ";
    const vector<string> &gp = m_particles.ghostParametersString();
    for (const string &param : gp) {
      cout << param << ", ";
    }
    cout << endl;
  }
#endif
  //--------------------------------------------------------------------------
  double nSec = timer.toc();
  if (isRoot)
    cout << "Time: " << nSec << "s" << endl;

  return 0;
}
//------------------------------------------------------------------------------
void PdSolver::solve() {
  if (isRoot)
    cout << "Starting solver" << endl;
  solver->solve();
}
//------------------------------------------------------------------------------
