#ifndef PD_PARTICLES_H
#define PD_PARTICLES_H

#include "particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
// Extending the particle struct for PD-data
//------------------------------------------------------------------------------

class PD_Particles: public Particles
{
protected:
    arma::mat m_r0;
    arma::mat m_v;
    arma::mat m_F;

    // For ADR
    arma::vec m_stableMass;
    arma::mat m_Fold;

    std::unordered_map<int, std::vector<std::pair<int, vector<double>>>> m_PdConnections;
    std::unordered_map<std::string, int> m_PdParameters;

public:
    PD_Particles();

    ~PD_Particles();


    virtual void initializeMatrices();

    void initializeADR();

    void setPdConnections(int id, vector<std::pair<int, std::vector<double>>> connections)
    {
        m_PdConnections[id] = connections;
    }

    std::vector<std::pair<int, std::vector<double>>> & pdConnections(int id)
    {
        return m_PdConnections.at(id);
    }

    arma::mat & r0()
    {
        return m_r0;
    }

    arma::mat & v()
    {
        return m_v;
    }

    arma::mat & F()
    {
        return m_F;
    }

    arma::vec & stableMass()
    {
        return m_stableMass;
    }

    arma::mat & Fold()
    {
        return m_Fold;
    }

    int getPdParamId(string paramId)
    {
        if(m_PdParameters.count(paramId) != 1)
        {
            cerr << "ERROR: accessing a PD_particles parameter that does not exist: "
                 << paramId << endl;
            throw ParameterDoesNotExist;
        }
        return m_PdParameters.at(paramId);
    }

    int registerPdParameter(string paramId, double value=0)
    {
        if(m_PdParameters.count(paramId) == 1)
        {
            cerr << "ERROR: PD-parameter already registerd: "
                 << paramId << endl;
            return m_PdParameters.at(paramId);
//            throw ParameterExist;
        }
        int pos = m_PdParameters.size();
        m_PdParameters[paramId] = pos;

        // Adding the new parameter to all connections
        for(int col=0;col<m_nParticles; col++)
        {
            for(auto &con:m_PdConnections[col])
            {
                con.second.push_back(value);
            }
        }

        return pos;
    }
};
//------------------------------------------------------------------------------
}
#endif // PD_PARTICLES_H
