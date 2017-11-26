#include "saveparticles.h"
#include "PDtools/Grid/grid.h"

#include <algorithm>
#include <cfloat>
#include "particles.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace PDtools
{
//------------------------------------------------------------------------------
// SaveParticles functions
//------------------------------------------------------------------------------
void SaveParticles::setGrid(Grid *mainGrid)
{
    m_mainGrid = mainGrid;
}
//------------------------------------------------------------------------------
void SaveParticles::writeToFile(Particles &particles, string savePath)
{
    initialize(particles);

    if(m_binary) {
        if(m_myRank == 0) {
            if(m_format == "xyz")
            {
                cerr << "Binary xyz not implemented" << endl;
                throw 10;
            }
            else if(m_format == "ply") {
                write_plyBinaryHeader(particles, savePath);
            }
            else if(m_format == "lmp") {
                write_lmpBinaryHeader(particles, savePath);
            }
        }
        else {
#if !USE_MPI
            savePath = savePath + to_string(m_myRank);
#endif
        }
        writeBinaryBody(particles, savePath);
    }
    else {
        if(m_myRank == 0) {
            if(m_format == "xyz") {
                write_xyzHeader(particles, savePath);
            }
            else if(m_format == "ply") {
                write_plyHeader(particles, savePath);
            }
            else if(m_format == "lmp") {
                write_lmpHeader(particles, savePath);
            }
        }
        else {
            savePath = savePath + to_string(m_myRank);
        }
        writeBody(particles, savePath);
    }
}
//------------------------------------------------------------------------------
void SaveParticles::setRankAndCores(int rank, int cores)
{
    m_myRank = rank;
    m_nCores = cores;
}
//------------------------------------------------------------------------------
void SaveParticles::initialize(Particles &particles)
{
    m_saveId = false;
    m_saveCoreId = false;
    m_saveCoordinates.clear();
    m_saveVelocities.clear();
    m_header.clear();

    // Checking for particle ids and positions
    for(const auto parameter_scale:m_saveParameters) {
        const string parameter = parameter_scale.first;
        const double scale = parameter_scale.second;

        if(parameter == "id") {
            m_saveId = true;
        }
        if(parameter == "coreId") {
            m_saveCoreId = true;
            m_header.push_back(parameter);
        }
        else if(parameter == "x") {
            m_saveCoordinates.push_back(pair<int, double>(0, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "y") {
            m_saveCoordinates.push_back(pair<int, double>(1, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "z") {
            m_saveCoordinates.push_back(pair<int, double>(2, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "v_x") {
            m_saveVelocities.push_back(pair<int, double>(0, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "v_y") {
            m_saveVelocities.push_back(pair<int, double>(1, scale));
            m_header.push_back(parameter);
        }
        else if(parameter == "v_z") {
            m_saveVelocities.push_back(pair<int, double>(2, scale));
            m_header.push_back(parameter);
        }
    }

    // Checking for other data
    m_dataParameters.clear();
    for(auto parameter_scale:m_saveParameters) {
        const string parameter = parameter_scale.first;
        const double scale = parameter_scale.second;

        if(particles.parameters().count(parameter)) {
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
    const ivec &colToId = particles.colToId();
    const int nParticles = particles.nParticles();

    ofstream outStream;
    outStream.open(savePath.c_str(), std::ofstream::out | std::ofstream::app);
    outStream.setf(ios::scientific);
    outStream.precision(14);

    for(int i=0; i<nParticles; i++) {
        const int id  = colToId(i);

        if(m_saveId) {
            outStream << id;
        }

        if(m_saveCoreId) {
            outStream << " " <<  m_myRank;
        }
        for(const auto & coord:m_saveCoordinates) {
            outStream << " " << particles.r()(i, coord.first)*coord.second;
        }

        for(const auto & coord:m_saveVelocities) {
            outStream << " " << particles.v()(i, coord.first)*coord.second;
        }

        for(const auto & parameter:m_dataParameters) {
            outStream << " " << particles.data()(i, parameter.first)*parameter.second;
        }

        outStream << endl;
    }
    //--------------------------------------------------------------------------
    outStream.close();

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);

    if(m_myRank <= 0) {
        // Writing the header
        ofstream of_fileConcatenated(savePath.c_str(), ios::out | ios::app );

        // Concatingating the files
        for(int node=1; node < m_nCores; node++) {
            string fName = savePath + to_string(node);
            std::ifstream if_node(fName, std::ios_base::binary);
            of_fileConcatenated << if_node.rdbuf();
            if_node.close();
            remove( fName.c_str() );
        }
        of_fileConcatenated.close();
    }
#endif
}
//--------------------------------------------------------------------------
void SaveParticles::writeBinaryBody(Particles &particles,
                                    const string &savePath)
{
    int nParticles = particles.nParticles();
#if USE_MPI
    string sPath = savePath;
    MPI_File binaryData;
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, &sPath[0], MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &binaryData);
    int nParticlesOnCPUs[m_nCores];
    MPI_Allgather(&nParticles, 1, MPI_INT, nParticlesOnCPUs, 1, MPI_INT, MPI_COMM_WORLD);
#else
    FILE* binaryData = fopen(savePath.c_str(), "a+");
#endif
    int headerOffset = 2*sizeof(long long int) + 7*sizeof(int) + 6*sizeof(double) + 3*sizeof(int);
    headerOffset /= sizeof(double);
    int bodyOffset = 0;
#if USE_MPI
    for(int i=0; i<m_myRank; i++)
    {
        bodyOffset += nParticlesOnCPUs[i];
    }
#endif

    const ivec &colToId = particles.colToId();
    const arma::mat & r = particles.r();
    const arma::mat & v = particles.v();
    const arma::mat & data = particles.data();
    int nColumns = 0;
    if(m_saveId) {
        nColumns++;
    }
    nColumns += m_header.size();
    int offset = headerOffset + bodyOffset*nColumns;


    for(int j=0; j<nParticles; j++) {
        const int id  = colToId(j);
        double buffer[nColumns];
        int i = 0;
        if(m_saveId) {
            buffer[i++] = id;
        }
        if(m_saveCoreId) {
            buffer[i++] = m_myRank;
        }

        for(const auto & coord:m_saveCoordinates) {
            const double r_i = r(j, coord.first)*coord.second;
            buffer[i++] = r_i;
        }

        for(const auto & coord:m_saveVelocities) {
            const double v_i = v(j, coord.first)*coord.second;
            buffer[i++] = v_i;
        }

        for(const auto & parameter:m_dataParameters) {
            const double dp = data(j, parameter.first)*parameter.second;
            buffer[i++] = dp;
        }

#if USE_MPI
        MPI_File_write_at(binaryData, offset*sizeof(double), &buffer, nColumns, MPI_DOUBLE, &status);
        offset += nColumns;
#else
        fwrite((char*)(&buffer), sizeof(double), nColumns, binaryData);
#endif
    }
    //--------------------------------------------------------------------------

#if USE_MPI
    MPI_File_close(&binaryData);
#else
    fclose(binaryData);

//    MPI_Barrier(MPI_COMM_WORLD);

//    if(m_myRank <= 0) {
//        // Writing the header
//        ofstream of_fileConcatenated(savePath.c_str(), ios::out | ios::app );

//        // Concatingating the files
//        for(int node=1; node < m_nCores; node++) {
//            string fName = savePath + to_string(node);
//            std::ifstream if_node(fName, std::ios_base::binary);
//            of_fileConcatenated << if_node.rdbuf();
//            if_node.close();
//            remove( fName.c_str() );
//        }
//        of_fileConcatenated.close();
//    }
#endif
}
//------------------------------------------------------------------------------
void SaveParticles::write_xyzHeader(const Particles &particles, const string &savePath)
{
    ofstream outStream;
    if(m_append) {
        outStream.open(savePath.c_str(), ios::app);
    }
    else {
        outStream.open(savePath.c_str());
    }
    outStream << particles.totParticles() << endl;

    // This program follows a convenction that the comments in a xyz-file
    // names the variables.

    outStream << "#";

    if(m_saveId) {
        outStream << " id";
    }

    for(string h:m_header) {
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
    if(m_append) {
        outStream.open(savePath.c_str(), ios::app);
    }
    else {
        outStream.open(savePath.c_str());
    }
    outStream << "ply" << endl;
    outStream << "format ascii 1.0" << endl;
    outStream << "element vertex " << particles.totParticles() << endl;

    if(m_saveId) {
        outStream << "property double id" << endl;
    }

    for(string h:m_header) {
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
    if(m_append) {
        outStream.open(savePath.c_str(), ios::app);
    }
    else {
        outStream.open(savePath.c_str());
    }

    const vector<pair<double, double>> & boundary = m_mainGrid->originalBoundary();
    const double scaling = m_saveCoordinates[0].second;
    const double xMin = scaling*boundary[0].first;
    const double xMax = scaling*boundary[0].second;
    const double yMin = scaling*boundary[1].first;
    const double yMax = scaling*boundary[1].second;
    const double zMin = scaling*boundary[2].first;
    const double zMax = scaling*boundary[2].second;


    outStream << "ITEM: TIMESTEP" << endl;
    outStream << m_timestep << endl;
    outStream << "ITEM: NUMBER OF ATOMS" << endl;
    outStream << particles.totParticles() << endl;
    outStream << "ITEM: BOX BOUNDS pp pp pp" << endl;
    outStream << xMin << ' ' << xMax << "\n";
    outStream << yMin << ' ' << yMax << "\n";
    outStream << zMin << ' ' << zMax << "\n";

    outStream << "ITEM: ATOMS";

    if(m_saveId) {
        outStream << " id";
    }

    for(string h:m_header) {
        outStream << " " << h;
    }
    outStream << endl;
}
//------------------------------------------------------------------------------
void SaveParticles::write_plyBinaryHeader(const Particles &particles,
                                          const string &savePath)
{
    ofstream outStream;
    if(m_append) {
        outStream.open(savePath.c_str(), ios::app);
    }
    else {
        outStream.open(savePath.c_str());
    }
    outStream << "ply" << endl;
    outStream << "format binary_little_endian 1.0" << endl;
    outStream << "element vertex " << particles.totParticles() << endl;

    if(m_saveId) {
        outStream << "property double id" << endl;
    }

    for(string h:m_header) {
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

    const vector<pair<double, double>> & boundary = m_mainGrid->originalBoundary();
    const double scaling = m_saveCoordinates[0].second;
    xMin = scaling*boundary[0].first;
    xMax = scaling*boundary[0].second;
    yMin = scaling*boundary[1].first;
    yMax = scaling*boundary[1].second;
    zMin = scaling*boundary[2].first;
    zMax = scaling*boundary[2].second;

    // Writing the header
    long long int currentTimeStep = m_timestep;
    long long int nParticles = particles.totParticles();
    int triclinic = 0.0;
    int nChunks = 1;

    int nColumns = 0;
    if(m_saveId) {
        nColumns++;
    }

    nColumns += m_header.size();

    int chunkLength = nParticles * nColumns;

    FILE* binaryData;
    if(m_append) {
        binaryData = fopen(savePath.c_str(), "a+");
    }
    else {
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
