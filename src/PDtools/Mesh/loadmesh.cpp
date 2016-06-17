#include "loadmesh.h"

#include <boost/algorithm/string.hpp>

//------------------------------------------------------------------------------
namespace PDtools
{
//------------------------------------------------------------------------------
PdMesh loadMesh2d(string path)
{
    // TODO: check path validity, format, etc
    PdMesh mesh = load_msh2d(path);
    return mesh;
}
//------------------------------------------------------------------------------
PdMesh load_msh2d(string loadPath)
{
    // Defintion of msh format:
    // http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php

    fstream data(loadPath, ios::in);
    string line;
    // Reading the number of particles
    getline(data, line);

    if(!boost::iequals(line, "$MeshFormat")) {

    }
    getline(data, line);

    std::vector<std::string> version_fileType_dataSize;
    boost::split(version_fileType_dataSize, line, boost::is_any_of(" "));

    const string meshVersion = version_fileType_dataSize[0];
    const string fileType = version_fileType_dataSize[1];
    const string dateSize = version_fileType_dataSize[2];

    bool inHeader = true;
    while(inHeader) {
        getline(data, line);
        inHeader = boost::iequals(line, "$EndMeshFormat") ? false:true;
    }

    bool inNodes = false;
    while(!inNodes) {
        getline(data, line);
        inNodes = boost::iequals(line, "$nodes");
    }
    getline(data, line);
    boost::trim(line);
    size_t nVertices = stoi(line);

    vector<int> verticeIds;
    vector<array<double, 3> > vertices;
    verticeIds.reserve(nVertices);
    vertices.reserve(nVertices);

    // Loading the ids and coordiantes
    for(size_t i=0; i<nVertices; i++) {
        getline(data, line);
        boost::trim_if(line, boost::is_any_of("\t "));
        std::vector<std::string> coordinateData;
        boost::split(coordinateData, line, boost::is_any_of("\t "), boost::token_compress_on);

        if(coordinateData.size() > 3) {
            verticeIds.push_back(stoi(coordinateData[0]));
            array<double, 3> coord = {stod(coordinateData[1]), stod(coordinateData[2]),  stod(coordinateData[3])};
            vertices.push_back(coord);
        } else {
            cerr << "Error reading mesh from file, vertices elements not matching." << endl;
            cerr << loadPath << endl;
        }
    }

    getline(data, line);
    if(!boost::iequals(line, "$EndNodes")) {
        cerr << "Error in reading mesh. Noes not ended correctly with $EndNodes" << endl;
    }

    bool inElements = false;
    while(!inElements) {
        getline(data, line);
        inElements = boost::iequals(line, "$Elements");
    }

    getline(data, line);
    boost::trim(line);
    size_t nElements = stoi(line);
    vector<vector<int>> elements;
    vector<int> elementIds;

    // Loading the connections
    for(size_t i=0; i<nElements; i++) {
        getline(data, line);
        boost::trim_if(line, boost::is_any_of("\t "));
        std::vector<std::string> elementData;
        boost::split(elementData, line, boost::is_any_of("\t "), boost::token_compress_on);

        const int elemendId_i = stoi(elementData[0]);
        const int elementType = stoi(elementData[1]);
        const int nTags = stoi(elementData[2]);
        vector<int> element;
        // Default tags: physicalEntity, geometricalEnetity

        switch(elementType) {
        case TRIANGLE:
            element =  {stoi(elementData[3 + nTags]),
                    stoi(elementData[3 + nTags + 1]),
                    stoi(elementData[3 + nTags + 2])};
        break;
        case QUADRANGLE:
            element = {stoi(elementData[3 + nTags]),
                       stoi(elementData[3 + nTags + 1]),
                       stoi(elementData[3 + nTags + 2]),
                       stoi(elementData[3 + nTags + 3])};
            break;
        default:
            cerr << "Error in reading element type.";
            break;
        }
        elementIds.push_back(elemendId_i);
        elements.push_back(element);
    }


    data.close();

    return PdMesh(vertices, verticeIds, elements, elementIds);
}
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
