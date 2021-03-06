#ifndef SAVEPARTICLES_H
#define SAVEPARTICLES_H

#include "config.h"

namespace PDtools {
class Particles;
class Grid;

//------------------------------------------------------------------------------
class SaveParticles {
private:
  string m_format;
  std::vector<std::pair<std::string, double>> m_saveParameters;
  bool m_binary = false;
  bool m_saveId;
  bool m_saveCoreId;
  vector<pair<int, double>> m_saveCoordinates;
  vector<pair<int, double>> m_saveVelocities;
  vector<pair<int, double>> m_dataParameters;
  vector<string> m_header;
  bool m_append = false;
  int m_timestep = 0;
  Grid *m_mainGrid;

  int m_myRank = 0;
  int m_nCores = 1;

public:
  SaveParticles();
  SaveParticles(string format, bool binary = false)
      : m_format(format), m_saveParameters({}), m_binary(binary) {
    ;
  }
  SaveParticles(string format,
                vector<pair<string, double>> saveParameters,
                bool binary = false)
      : m_format(format), m_saveParameters(saveParameters), m_binary(binary) {
    ;
  }

  void writeToFile(Particles &particles, string savePath);
  void append(bool append) { m_append = append; }
  void binary(bool binary) { m_binary = binary; }
  void setTimestep(int timestep) { m_timestep = timestep; }
  void setRankAndCores(int rank, int cores);
  void setGrid(Grid *mainGrid);

private:
  void initialize(Particles &particles);

  void writeBody(Particles &particles, const string &savePath);
  void writeBinaryBody(Particles &particles, const string &savePath);
  void write_xyzHeader(const Particles &particles, const string &savePath);
  void write_plyHeader(const Particles &particles, const string &savePath);
  void write_lmpHeader(const Particles &particles, const string &savePath);
  void write_plyBinaryHeader(const Particles &particles,
                             const string &savePath);
  void write_lmpBinaryHeader(Particles &particles, const string &savePath);
};
}
//------------------------------------------------------------------------------
#endif // SAVEPARTICLES_H
