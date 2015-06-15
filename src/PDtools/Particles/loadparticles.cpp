#include "loadparticles.h"

#include "pd_particles.h"
#include "particles.h"

namespace PDtools
{
template class T_LoadParticles<Particles>;
template class T_LoadParticles<PD_Particles>;

//------------------------------------------------------------------------------
template <class T_particles>T_LoadParticles<T_particles>::T_LoadParticles()
{
}
//------------------------------------------------------------------------------
template <class T_particles>
T_particles T_LoadParticles<T_particles>::load(string loadPath,
                               string format,
                               bool bin,
                               unordered_map<string, int> loadParameters)
{
    m_format = format;
    m_binary = bin;
    T_particles particles;

    if(m_binary)
    {
        FILE* binaryData = fopen(loadPath.c_str(), "rb");

        if (!binaryData)
        {
          cerr << "ERROR: Could not open " << loadPath << endl;
        }

        if(m_format == "xyz")
        {
            cerr << "Binary xyz not implemented" << endl;
            throw 10;
        }
        else if(m_format == "ply")
        {
            loadParameters = read_plyBinaryHeader(binaryData, loadParameters);
        }
        else if(m_format == "lmp")
        {
            read_lmpBinaryHeader(binaryData, loadParameters);
        }
        else
        {
            cerr << "Format: '" << m_format << "' not supported" << endl;
            throw 10;
        }

        loadBinaryBody(particles, binaryData, loadParameters);
        fclose(binaryData);
    }
    else
    {
        fstream data(loadPath, ios::in);
        unordered_map <std::string, int> parameters;

        if(m_format == "xyz")
        {
            parameters = read_xyzHeader(data);
        }
        else if(m_format == "ply")
        {
            cerr << "ply not implemented" << endl;
            throw 10;
        }
        else if(m_format == "lmp")
        {
            cerr << "lmp not implemented" << endl;
            throw 10;
        }
        else
        {
            cerr << "Format: '" << m_format << "' not supported" << endl;
            throw 10;
        }

        loadBody(particles, data, parameters);
        data.close();
    }

    return particles;
}
//------------------------------------------------------------------------------
template <class T_particles>
unordered_map <std::string, int> T_LoadParticles<T_particles>::read_xyzHeader(fstream &data)
{
    string line;

    // Reading the number of particles
    getline(data, line);

    m_nParticles = stoi(line);

    // Reading the comments
    // This program follows a convenction that the comments
    // in a xyz-file must name the variables.
    getline(data, line);
    regex rr("([A-Za-z_]+)");
    sregex_iterator next(line.begin(), line.end(), rr);
    sregex_iterator end;

    unordered_map <std::string, int> parameters;
    int position = 0;
    while(next != end)
    {
        smatch match = *next;
        parameters[match.str()] = position;
        position++;
        next++;
    }

    return parameters;
}
//------------------------------------------------------------------------------
template <class T_particles>
unordered_map<string, int> T_LoadParticles<T_particles>::read_lmpBinaryHeader(
        FILE *data,
        unordered_map<string, int> parameters)
{
    (void) parameters;
    int timeStep;
    int nParticles;
    int triclinic = false;
    int nChunks;
    int nColumns;
    int chunkLength;
    int boundary[3][2];
    double boundaryCoordinates[3][2];
    double shear[3];

    if(m_useLegacyFormat)
    {
        fread(&timeStep, sizeof(int), 1, data);
        fread(&nParticles,  sizeof(int), 1, data);
    }
    else
    {
        long long int tmp_currentTimeStep, tmp_nParticles;
        fread(&tmp_currentTimeStep, sizeof(long long int), 1, data);
        fread(&tmp_nParticles,  sizeof(long long int), 1, data);
        timeStep = tmp_currentTimeStep;
        nParticles = tmp_nParticles;
        fread(&triclinic,   sizeof(int), 1, data);
        fread(&boundary[0][0], 6*sizeof(int), 1, data);
    }
    fread(&boundaryCoordinates[0][0], 6*sizeof(double), 1, data);

    if(m_useLegacyFormat)
    {
        fread(&shear[0], 3*sizeof(double), 1, data);
    }

    if (triclinic)
    {
        cerr << "ERROR: triclinic not supported in LAMMPS binary file" << endl;
        int dump[3];
        fread(&dump[0], 3*sizeof(int), 1 ,data);
    }


    fread(&nColumns,    sizeof(int), 1, data);
    fread(&nChunks,     sizeof(int), 1, data);
    fread(&chunkLength, sizeof(int), 1, data);

    m_timeStep = timeStep;
    m_nParticles = nParticles;
    m_nColumns = nColumns;
    m_chunkLength = chunkLength;

    unordered_map<string, int> empty;
    return empty;
}
//------------------------------------------------------------------------------
template <class T_particles>
unordered_map<string, int> T_LoadParticles<T_particles>::read_plyBinaryHeader(
        FILE *data,
        unordered_map<string, int> parameters)
{
    parameters = {};
    string line = readBinaryLine(data);

    if(line != "ply")
    {
        cerr << "Error: loading corrupt ply-file:";
        throw 20;
    }

    line = readBinaryLine(data);
    if(line != "format binary_little_endian 1.0")
    {
        cerr << "Error: error loading binary ply-file:";
        throw 20;
    }

    int position = 0;
    int nParticles = 0;
    while(line != "end_header")
    {
        line = readBinaryLine(data);
        vector<string> lineSplit;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::split(lineSplit, line, boost::is_any_of("\t "), boost::token_compress_on);

        if(lineSplit[0] == "comment")
            continue;

        if(lineSplit[0] == "element")
        {
            if(lineSplit[1] == "vertex")
            {
                nParticles = stoi(lineSplit[2]);
            }
            else
            {
                cerr << "Error: only vertex element supported in ply file" << endl;
                cerr << "found: " << lineSplit[1] << endl;
                throw 20;
            }
        }

        if(lineSplit[0] == "property")
        {
            if(lineSplit[1] == "double")
            {
                parameters[lineSplit[2]] = position;
                position++;
            }
            else
            {
                cerr << "Error: in ply 'property'' only 'double' is supported" << endl;
                cerr << "found: " << lineSplit[1] << endl;
                throw 20;
            }
        }
    }

    m_timeStep = 0;
    m_nColumns = parameters.size();
    m_nParticles = nParticles;
    m_chunkLength = nParticles*parameters.size();

    return parameters;
}
//------------------------------------------------------------------------------
template <class T_particles>
void T_LoadParticles<T_particles>::loadBody(T_particles &particles, fstream &rawData,
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
    //--------------------------------------------------------------------------
    // Creating the data matrix
    particles.initializeMatrices();
    unordered_map<int, int> & pIds = particles.pIds();
    arma::ivec & get_id = particles.get_id();
    arma::mat & r = particles.r();
    arma::mat & data = particles.data();

    // Reading all the data from file
    for(int i=0; i<particles.nParticles(); i++)
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
            r(pc.first, i) = stod(lineSplit[pc.second]);
        }
        for(pair<int, int> dfc:data_config_mapping)
        {
            data(i, dfc.first) = stod(lineSplit[dfc.second]);
        }
    }
}
//------------------------------------------------------------------------------
template <class T_particles>
void T_LoadParticles<T_particles>::loadBinaryBody(T_particles &particles,
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
    if(parameters.count("id") > 0)   {
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

    //--------------------------------------------------------------------------
    // Creating the data matrix
    particles.initializeMatrices();

    int nColumns = m_nColumns;
    unordered_map<int, int> & pIds = particles.pIds();
    arma::ivec & get_id = particles.get_id();
    arma::mat & r = particles.r();
    arma::mat & data = particles.data();

    // Reading all the data from file
    for(int i=0; i<particles.nParticles(); i++)
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
            r(pc.first, i) = line[pc.second];
        }
        for(pair<int, int> dfc:data_config_mapping)
        {
            data(i, dfc.first) = line[dfc.second];
        }
    }
}
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
