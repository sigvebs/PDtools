#include "saveparticles.h"

#include <algorithm>
#include <cfloat>
#include "particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
// SaveParticles functions
//------------------------------------------------------------------------------
void SaveParticles::writeToFile(Particles &particles, const string &savePath)
{
    initialize(particles);

    if(m_binary)
    {
        if(m_format == "xyz")
        {
            cerr << "Binary xyz not implemented" << endl;
            throw 10;
        }
        else if(m_format == "ply")
        {
            write_plyBinaryHeader(particles, savePath);
        }
        else if(m_format == "lmp")
        {
            write_lmpBinaryHeader(particles, savePath);
        }

        writeBinaryBody(particles, savePath);
    }
    else
    {

        if(m_format == "xyz")
        {
            write_xyzHeader(particles, savePath);
        }
        else if(m_format == "ply")
        {
            write_plyHeader(particles, savePath);
        }
        else if(m_format == "lmp")
        {
            write_lmpHeader(particles, savePath);
        }

        writeBody(particles, savePath);
    }
}
//------------------------------------------------------------------------------
void SaveParticles::initialize(Particles &particles)
{
    // The default is that all data-parameters are saved
    /*
    if(m_saveParameters.empty())
    {
        m_saveParameters.push_back("id");
        m_saveParameters.push_back("x");
        m_saveParameters.push_back("y");
        m_saveParameters.push_back("z");

        for(auto id_pos:particles.parameters())
        {
            m_saveParameters.push_back(id_pos.first);
        }
    }*/

    m_saveId = false;
    m_saveCoordinates.clear();
    m_header.clear();

    // Checking for particle ids and positions
    for(const auto parameter_scale:m_saveParameters)
    {
        const string parameter = parameter_scale.first;
        const double scale = parameter_scale.second;

        if(parameter == "id")
        {
            m_saveId = true;
        }
        else if(parameter == "x")
        {
            m_saveCoordinates.push_back(pair<int, double>(0, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "y")
        {
            m_saveCoordinates.push_back(pair<int, double>(1, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "z")
        {
            m_saveCoordinates.push_back(pair<int, double>(2, scale));
            m_header.push_back(parameter);
        }
    }

    // Checking for other data
    m_dataParameters.clear();
    for(auto parameter_scale:m_saveParameters)
    {
        const string parameter = parameter_scale.first;
        const double scale = parameter_scale.second;

        if(particles.parameters().count(parameter))
        {
            std::pair<int, double> p(particles.parameters().at(parameter), scale);
            m_dataParameters.push_back(p);
            m_header.push_back(parameter);
        }
    }
}
//------------------------------------------------------------------------------
void SaveParticles::writeBody(Particles &particles,
                             const string &savePath)
{
    ofstream outStream;
    outStream.open(savePath.c_str(), std::ofstream::out | std::ofstream::app);
    outStream.setf(ios::scientific);
    outStream.precision(14);

    for(auto id_pos:particles.pIds())
    {
        int id  = id_pos.first + 1;
        int j = id_pos.second;

        if(m_saveId)
        {
            outStream << id;
        }

        for(const auto & coord:m_saveCoordinates)
        {
            outStream << " " << particles.r()(j, coord.first)*coord.second;
        }

        for(const auto & parameter:m_dataParameters)
        {
            outStream << " " << particles.data()(j, parameter.first)*parameter.second;
        }

        outStream << endl;
    }
    //--------------------------------------------------------------------------
    outStream.close();
}
//--------------------------------------------------------------------------
void SaveParticles::writeBinaryBody(Particles &particles,
                                   const string &savePath)
{
    FILE* binaryData = fopen(savePath.c_str(), "a+");
    const arma::mat & r = particles.r();
    const arma::mat & data = particles.data();

    for(auto id_pos:particles.pIds())
    {
        const double id  = id_pos.first + 1;
        int j = id_pos.second;

        if(m_saveId)
        {
            fwrite((char*)(&id), sizeof(double), 1, binaryData);
        }

        for(const auto & coord:m_saveCoordinates)
        {
            const double r_i = r(j, coord.first)*coord.second;
            fwrite((char*)(&r_i), sizeof(double), 1, binaryData);
        }

        for(const auto & parameter:m_dataParameters)
        {
            const double dp = data(j, parameter.first)*parameter.second;
            fwrite((char*)(&dp), sizeof(double), 1, binaryData);
        }
    }
    //--------------------------------------------------------------------------
    fclose(binaryData);
}
//------------------------------------------------------------------------------
void SaveParticles::write_xyzHeader(const Particles &particles, const string &savePath)
{
    ofstream outStream;
    if(m_append)
    {
        outStream.open(savePath.c_str(), ios::app);
    }
    else
    {
        outStream.open(savePath.c_str());
    }
    outStream << particles.nParticles() << endl;

    // This program follows a convenction that the comments in a xyz-file
    // names the variables.

    outStream << "#";

    if(m_saveId)
    {
        outStream << " id";
    }

    for(string h:m_header)
    {
        outStream << " " << h;
    }
    outStream << endl;
    outStream.close();
}
//------------------------------------------------------------------------------
void SaveParticles::write_plyHeader(const Particles &particles,
                                   const string &savePath)
{
    ofstream outStream;
    if(m_append)
    {
        outStream.open(savePath.c_str(), ios::app);
    }
    else
    {
        outStream.open(savePath.c_str());
    }
    outStream << "ply" << endl;
    outStream << "format ascii 1.0" << endl;
    outStream << "element vertex " << particles.nParticles() << endl;

    if(m_saveId)
    {
        outStream << "property double id" << endl;
    }

    for(string h:m_header)
    {
        outStream << "property float " << h << endl;
    }
    outStream << "end_header" << endl;
    outStream.close();
}
//------------------------------------------------------------------------------
void SaveParticles::write_lmpHeader(const Particles &particles,
                                   const string &savePath)
{
    ofstream outStream;
    if(m_append)
    {
        outStream.open(savePath.c_str(), ios::app);
    }
    else
    {
        outStream.open(savePath.c_str());
    }
    outStream << "ITEM: TIMESTEP" << endl;
    outStream << m_timestep << endl;
    outStream << "ITEM: NUMBER OF ATOMS" << endl;
    outStream << particles.nParticles() << endl;

    outStream << "ITEM: ATOMS";

    if(m_saveId)
    {
        outStream << " id";
    }

    for(string h:m_header)
    {
        outStream << " " << h;
    }
    outStream << endl;
}
//------------------------------------------------------------------------------
void SaveParticles::write_plyBinaryHeader(const Particles &particles,
                                         const string &savePath)
{
    ofstream outStream;
    if(m_append)
    {
        outStream.open(savePath.c_str(), ios::app);
    }
    else
    {
        outStream.open(savePath.c_str());
    }
    outStream << "ply" << endl;
    outStream << "format binary_little_endian 1.0" << endl;
    outStream << "element vertex " << particles.nParticles() << endl;

    if(m_saveId)
    {
        outStream << "property double id" << endl;
    }

    for(string h:m_header)
    {
        outStream << "property double " << h << endl;
    }
    outStream << "end_header" << endl;
    outStream.close();
}
//------------------------------------------------------------------------------
void SaveParticles::write_lmpBinaryHeader(Particles &particles,
                                         const string &savePath)
{
    // Finding the domain boundaries
    double xMin = DBL_MAX;
    double xMax = -DBL_MAX;
    double yMin = DBL_MAX;
    double yMax = -DBL_MAX;
    double zMin = DBL_MAX;
    double zMax = -DBL_MAX;

    int b_xMin = 2;
    int b_xMax = 2;
    int b_yMin = 2;
    int b_yMax = 2;
    int b_zMin = 2;
    int b_zMax = 2;

    for(auto id_pos:particles.pIds())
    {
        const int j = id_pos.second;
        xMin = xMin < particles.r()(j, 0) ? xMin : particles.r()(j, 0);
        xMax = xMax > particles.r()(j, 0) ? xMax : particles.r()(j, 0);
        yMin = yMin < particles.r()(j, 1) ? yMin : particles.r()(j, 1);
        yMax = yMax > particles.r()(j, 1) ? yMax : particles.r()(j, 1);
        zMin = zMin < particles.r()(j, 2) ? zMin : particles.r()(j, 2);
        zMax = zMax > particles.r()(j, 2) ? zMax : particles.r()(j, 2);
    }

    // Writing the header
    long long int currentTimeStep = m_timestep;
    long long int nParticles = particles.nParticles();
    int triclinic = 0.0;
    int nChunks = 1;

    int nColumns = 0;
    if(m_saveId)
    {
        nColumns++;
    }

    nColumns += m_header.size();

    int chunkLength = nParticles * nColumns;

    FILE* binaryData;
    if(m_append)
    {
        binaryData = fopen(savePath.c_str(), "a+");
    }
    else
    {
        binaryData = fopen(savePath.c_str(), "wb");
    }

    fwrite(&currentTimeStep, sizeof(long long int), 1, binaryData);
    fwrite(&nParticles, sizeof(long long int), 1, binaryData);
    fwrite(&triclinic, sizeof(int), 1, binaryData);
    fwrite(&b_xMin, sizeof(int), 1, binaryData);
    fwrite(&b_xMax, sizeof(int), 1, binaryData);
    fwrite(&b_yMin, sizeof(int), 1, binaryData);
    fwrite(&b_yMax, sizeof(int), 1, binaryData);
    fwrite(&b_zMin, sizeof(int), 1, binaryData);
    fwrite(&b_zMax, sizeof(int), 1, binaryData);
    fwrite(&xMin, sizeof(double), 1, binaryData);
    fwrite(&xMax, sizeof(double), 1, binaryData);
    fwrite(&yMin, sizeof(double), 1, binaryData);
    fwrite(&yMax, sizeof(double), 1, binaryData);
    fwrite(&zMin, sizeof(double), 1, binaryData);
    fwrite(&zMax, sizeof(double), 1, binaryData);
    fwrite(&nColumns, sizeof(int), 1, binaryData);
    fwrite(&nChunks, sizeof(int), 1, binaryData);
    fwrite(&chunkLength, sizeof(int), 1, binaryData);

    fclose(binaryData);
}
//------------------------------------------------------------------------------
}
