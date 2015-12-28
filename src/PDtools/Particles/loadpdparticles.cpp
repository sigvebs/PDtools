#include "loadpdparticles.h"

#include "pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
LoadPdParticles::LoadPdParticles()
{
}
//------------------------------------------------------------------------------
void LoadPdParticles::loadBody(PD_Particles &particles,
                             fstream &rawData,
                             unordered_map<string, int> parameters)
{
    particles.nParticles(m_nParticles);
    string line;

    //--------------------------------------------------------------------------
    // Storing only non-basic parameters in the parameters
    int counter = 0;
    vector<pair<int, int>> data_config_mapping;
    for(auto param:parameters)
    {
        bool found = false;
        for(string basic_parameter:m_PdBasicParameters)
        {
            if(param.first == basic_parameter)
            {
                found = true;
            }
        }

        if(!found)
        {
            particles.parameters()[param.first] = counter;
            data_config_mapping.push_back(pair<int,int>(counter, param.second));
            counter++;
        }
    }

    //--------------------------------------------------------------------------
    bool idIsset = false;
    int idPos = 0;
    if(parameters.count("id") > 0)
    {
        idIsset = true;
        idPos = parameters["id"];
    }

    vector<pair<int, int>> position_config;
    int dim = 0;
    if(parameters.count("x") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(0, parameters["x"]));
    }
    if(parameters.count("y") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(1, parameters["y"]));
    }
    if(parameters.count("z") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(2, parameters["z"]));
    }
    particles.dim(dim);

    vector<pair<int, int>> velocity_config;
    if(parameters.count("v_x") > 0)
    {
        velocity_config.push_back(pair<int, int>(0, parameters["v_x"]));
    }
    if(parameters.count("v_y") > 0)
    {
        velocity_config.push_back(pair<int, int>(1, parameters["v_y"]));
    }
    if(parameters.count("v_z") > 0)
    {
        velocity_config.push_back(pair<int, int>(2, parameters["v_z"]));
    }

    //--------------------------------------------------------------------------
    // Creating the data matrix
    particles.initializeMatrices();

    unordered_map<int, int> & pIds = particles.pIds();
    arma::ivec & get_id = particles.get_id();
    arma::mat & r = particles.r();
    arma::mat & v = particles.v();
    arma::mat & data = particles.data();


    // Reading all the data from file
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        vector<string> lineSplit;
        getline(rawData, line);
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::split(lineSplit, line, boost::is_any_of("\t "), boost::token_compress_on);

        // Collecting the data
        if(idIsset)
        {
            pIds[stoi(lineSplit[idPos])] = i;
            get_id[i] = stoi(lineSplit[idPos]);
        }
        else
        {
            pIds[i] = i;
            get_id[i] = i;
        }

        for(pair<int, int> pc:position_config)
        {
            r(i, pc.first) = stod(lineSplit[pc.second]);
        }
        for(pair<int, int> vc:velocity_config)
        {
            v(i, vc.first) = stod(lineSplit[vc.second]);
        }
        for(pair<int, int> dfc:data_config_mapping)
        {
            data(i, dfc.first) = stod(lineSplit[dfc.second]);
        }
    }
}
//------------------------------------------------------------------------------
void LoadPdParticles::loadBinaryBody(PD_Particles &particles,
                                     FILE *rawData,
                                     unordered_map<string, int> parameters)
{
    particles.nParticles(m_nParticles);
    // Storing only non-basic parameters in the parameters
    int counter = 0;
    vector<pair<int, int>> data_config_mapping;
    for(auto param:parameters)
    {
        bool found = false;
        for(string basic_parameter:basicParameters)
        {
            if(param.first == basic_parameter)
            {
                found = true;
            }
        }

        if(!found)
        {
            particles.parameters()[param.first] = counter;
            data_config_mapping.push_back(pair<int,int>(counter, param.second));
            counter++;
        }
    }

    //--------------------------------------------------------------------------
    bool idIsset = false;
    int idPos = 0;
    if(parameters.count("id") > 0)
    {
        idIsset = true;
        idPos = parameters["id"];
    }

    vector<pair<int, int>> position_config;
    int dim = 0;
    if(parameters.count("x") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(0, parameters["x"]));
    }
    if(parameters.count("y") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(1, parameters["y"]));
    }
    if(parameters.count("z") > 0)
    {
        dim++;
        position_config.push_back(pair<int, int>(2, parameters["z"]));
    }
    particles.dim(dim);

    vector<pair<int, int>> velocity_config;
    if(parameters.count("v_x") > 0)
    {
        velocity_config.push_back(pair<int, int>(0, parameters["v_x"]));
    }
    if(parameters.count("v_y") > 0)
    {
        velocity_config.push_back(pair<int, int>(1, parameters["v_y"]));
    }
    if(parameters.count("v_z") > 0)
    {
        velocity_config.push_back(pair<int, int>(2, parameters["v_z"]));
    }

    //--------------------------------------------------------------------------
    // Creating the data matrix
    particles.initializeMatrices();

    unordered_map<int, int> & pIds = particles.pIds();
    arma::ivec & get_id = particles.get_id();
    arma::mat & r = particles.r();
    arma::mat & v = particles.v();
    arma::mat & data = particles.data();

    int nColumns = m_nColumns;

    // Reading all the data from file
    for(unsigned int i=0; i<particles.nParticles(); i++)
    {
        double line[nColumns];
        fread(&line[0], nColumns*sizeof(double), 1, rawData);

        // Collecting the data
        if(idIsset)
        {
            pIds[int(line[idPos])] = i;
            get_id[i] = int(line[idPos]);
        }
        else
        {
            pIds[i] = i;
            get_id[i] = i;
        }

        for(pair<int, int> pc:position_config)
        {
            r(i, pc.first) = line[pc.second];
        }
        for(pair<int, int> vc:velocity_config)
        {
            v(i, vc.first) = line[vc.second];
        }
        for(pair<int, int> dfc:data_config_mapping)
        {
            data(i, dfc.first) = line[dfc.second];
        }
    }
}
//------------------------------------------------------------------------------
PD_Particles load_pd(string loadPath)
{
    LoadPdParticles loadParticles;
    vector<string> lineSplit;

    boost::split(lineSplit, loadPath, boost::is_any_of("."), boost::token_compress_on);
    const string type = lineSplit.back();

    PD_Particles particles = loadParticles.load(loadPath, type);

    particles.type(type);
    return particles;
}
//------------------------------------------------------------------------------
PD_Particles load_pd(string loadPath, unordered_map<string, double> particleParameters)
{
    PD_Particles particles = load_pd(loadPath);

    if(!particleParameters.empty())
    {
        for(auto pm:particleParameters)
        {
            particles.registerParameter(pm.first, pm.second);
        }
    }
    return particles;
}
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------

