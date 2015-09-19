#ifndef SAVEPARTICLES_H
#define SAVEPARTICLES_H

#include <vector>
#include <fstream>
using namespace std;

namespace PDtools
{
// Forward declerations
class Particles;

//------------------------------------------------------------------------------
class SaveParticles
{
private:
    string m_format;
    vector<string> m_saveParameters;
    bool m_binary = false;
    bool m_saveId;
    vector<int> m_saveCoordinates;
    vector<int> m_dataParameters;
    vector<string> m_header;
    bool m_append = false;
    int m_timestep = 0;

public:
    SaveParticles();
    SaveParticles(string format,
                 bool binary = false):
        m_format(format),
        m_saveParameters({}),
        m_binary(binary)
    {
        ;
    }
    SaveParticles(string format,
                 vector<string> saveParameters,
                 bool binary = false):
        m_format(format),
        m_saveParameters(saveParameters),
        m_binary(binary)
    {
        ;
    }
    void writeToFile(Particles &particles, const string &savePath);

    void append(bool append)
    {
        m_append = append;
    }
    void binary(bool binary)
    {
        m_binary = binary;
    }

    void setTimestep(int timestep)
    {
        m_timestep = timestep;
    }

private:
    void initialize(Particles &particles);

    void writeBody(Particles &particles,
                   const string &savePath);
    void writeBinaryBody(Particles &particles,
                         const string &savePath);

    void write_xyzHeader(const Particles &particles,
                         const string &savePath);
    void write_plyHeader(const Particles &particles,
                         const string &savePath);
    void write_lmpHeader(const Particles &particles,
                         const string &savePath);
    void write_plyBinaryHeader(const Particles &particles,
                         const string &savePath);
    void write_lmpBinaryHeader(Particles &particles,
                         const string &savePath);
};
}
//------------------------------------------------------------------------------
#endif // SAVEPARTICLES_H
