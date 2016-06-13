#ifndef LOADPARTICLES_H
#define LOADPARTICLES_H

#include <armadillo>
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace PDtools
{
// Forward declerations
class Particles;
class Grid;

//------------------------------------------------------------------------------
class LoadParticles
{
protected:
    Grid * m_grid = nullptr;
    int m_nParticles;
    string m_format;
    bool m_binary;
    const vector<string> basicParameters = {"id", "x", "y", "z"};
    bool m_useLegacyFormat = false;
    vector<pair<int, int>> data_config_mapping;
    int m_timeStep = 0;
    int m_nColumns;
    int m_chunkLength;

    virtual unordered_map <string, int>
    read_xyzHeader(fstream &data);

    virtual unordered_map <string, int>
    read_lmpHeader(fstream &data);

    virtual unordered_map <string, int>
    read_lmpBinaryHeader(FILE *data, unordered_map<string, int> parameters);

    virtual unordered_map <string, int>
    read_plyBinaryHeader(FILE *data, unordered_map<string, int> parameters);

    virtual void
    loadBody(Particles &particles, fstream &data, unordered_map <string, int> parameters);

    virtual void
    loadBinaryBody(Particles &particles, FILE *data, unordered_map <string, int> parameters);

public:
    LoadParticles();
    ~LoadParticles();

    Particles
    load(string loadPath,
         string format,
         bool bin=false,
         unordered_map<string, int> loadParameters={{}});

    void
    setGrid(Grid & grid)
    {
        m_grid = &grid;
    }

    void
    useLegacyFormat(bool ulf)
    {
        m_useLegacyFormat = ulf;
    }
};

//------------------------------------------------------------------------------
inline string readBinaryLine(FILE *data)
{
    string line = "";
    char inChar;
    int counter = 0;

    fread(&inChar, sizeof(char), 1, data);

    while(inChar != '\n' || counter >= 10000)
    {
        line += inChar;
        fread(&inChar, sizeof(char), 1, data);
    }

    if(counter >= 10000)
    {
        cerr << "Error reading binary line" << endl;
        throw 20;
    }

    return line;
}
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
#endif // LOADPARTICLES_H

//    // Regex to find a number on scientific form
//    regex re("((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?((e|E)((\\+|-)?)[[:digit:]]+)?");

// Example of xyz-file
//10
//# id x y v_x v_y volume
//0	0.313150    0.597834	0.165347	0.265347	8.76566e-06
//1	0.932484	0.519061	0.001397   	0.101397   	8.22899e-06
//2	0.731905	0.547958	0.212907	0.112907	8.0501e-06
//3	0.650592	0.343276	0.529744	0.329744	8.31843e-06
//4	0.109496	0.449608	0.195344	0.695344	7.78176e-06
//5	0.532164	0.704623	0.504753	0.0504753	6.44008e-06
//6	0.685139	0.468903	0.236449	0.536449	7.33453e-06
//7	0.678746	0.390759	0.546610	0.346610	8.0501e-06
//8	0.384222	0.581632	0.681405	0.168145	8.31843e-06
//9	0.947418	0.401296	0.665064	0.166564	7.0662e-06
//10	0.294880	0.368566	0.110014	0.111014	8.13954e-06
