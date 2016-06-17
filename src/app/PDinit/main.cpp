#include <iostream>
#include <armadillo>
#include <memory>
#include <libconfig.h++>

#include "PDtools.h"
#include "Mesh/pdmesh.h"
#include "Mesh/loadmesh.h"
#include "Mesh/meshtopdpartices.h"
#include "PdFunctions/pdfunctions.h"

using namespace std;
using namespace PDtools;
using namespace libconfig;

int main(int argc, char** argv)
{

    (void) argv;
    if(argc != 2)
    {
        cerr << "usage: peridyn [path to config file]" << endl;
        return 1;
    }
#ifdef USE_MPI
    MPI::Init (argc, argv);
#endif
    const int dim = 2;
    int m_myRank = 0;
    int m_nCores = 1;

    string cfgPath = argv[1];
    libconfig::Config m_cfg;
    m_cfg.readFile(cfgPath.c_str());

    double L0 = 1.;
    double dxdydz[M_DIM];
    Setting &cfg_domain = m_cfg.lookup("domain");
    vector<pair<double,double>> domain;
    for(int d=0; d<M_DIM; d++)
    {
        pair<double, double> bound(cfg_domain[2*d], cfg_domain[2*d+1]);
        bound.first /= L0;
        bound.second /= L0;
        dxdydz[d] = bound.second - bound.first;
        domain.push_back(bound);
    }

    // Setting periodicity
    arma::ivec3 periodicBoundaries = {0, 0, 0};

    try
    {
        Setting &cfg_periodic= m_cfg.lookup("periodic");
        for(int d=0; d<M_DIM; d++)
        {
            periodicBoundaries(d) = (int) cfg_periodic[d];
        }
    }
    catch(libconfig::SettingNotFoundException s)
    {
        if(m_myRank == 0)
            cout << "No periodicity set." << endl;
    }
    //--------------------------------------------------------------------------
    // Setting the grid
    //--------------------------------------------------------------------------
    Grid m_grid;

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

    double gridspacing = 1.25*(delta + 0.5*lc);
    m_grid = Grid(domain, gridspacing, periodicBoundaries);
    m_grid.setIdAndCores(m_myRank, m_nCores);
    m_grid.dim(dim);
    m_grid.initialize();
    m_grid.setMyGridpoints();
    m_grid.setInitialPositionScaling(L0);
    if(m_myRank == 0)
    {
        vector<int> cpuConfig = m_grid.nCpuGrid();
        cout << "Running core configuration: ";
        for(int n:cpuConfig)
            cout << n << " ";
        cout << endl;
    }
    //--------------------------------------------------------------------------
    // Loading mesh/particle configuration
    //--------------------------------------------------------------------------
    string meshPath = (const char *) m_cfg.lookup("particlesPath");
    cout << meshPath << endl;
    PD_Particles particles;

    string fileType = getFileEnding(meshPath);
    if(boost::iequals(fileType, "msh")) {
        int quadratureDegree = 1;
        m_cfg.lookupValue("quadratureDegree", quadratureDegree);
        PdMesh msh = loadMesh2d(meshPath);
        particles = convertMshToPdParticles(dim, quadratureDegree, msh, m_grid);
    }


    //--------------------------------------------------------------------------
    // Saving the data
    //--------------------------------------------------------------------------
    vector<pair<string, double>> saveparam_scale;
    saveparam_scale.push_back(pair<std::string, double>("id", 1.));
    if(dim >= 1)
        saveparam_scale.push_back(pair<std::string, double>("x", 1.));
    if(dim >= 2)
        saveparam_scale.push_back(pair<std::string, double>("y", 1.));
    if(dim >= 3)
        saveparam_scale.push_back(pair<std::string, double>("z", 1.));

    saveparam_scale.push_back(pair<std::string, double>("volume", 1.));

    string saveParticlesPath = "geometry.xyz";
    SaveParticles *saveParticles = new SaveParticles("xyz", saveparam_scale, false);
    saveParticles->writeToFile(particles, saveParticlesPath);


    cout << "Done initPeridyn." << endl;
#ifdef USE_MPI
    MPI::Finalize( );
#endif
    return 0;
}
