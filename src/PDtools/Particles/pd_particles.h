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
    mat m_r0;
    mat m_v;
    mat m_F;
    mat m_b;
    mat m_u;

    // For ADR
    vec m_stableMass;
    mat m_Fold;

    unordered_map<int, vector<pair<int, vector<double>>>> m_PdConnections;
    unordered_map<string, int> m_PdParameters;

public:
    PD_Particles();

    ~PD_Particles();

    virtual void
    initializeMatrices();

    void
    initializeADR();

    void
    initializeBodyForces();

    void
    setPdConnections(int id, vector<pair<int, vector<double>>> connections)
    {
        m_PdConnections[id] = connections;
    }

    vector<pair<int, vector<double>>> &
    pdConnections(int id)
    {
        return m_PdConnections.at(id);
    }

    mat &
    r0()
    {
        return m_r0;
    }

    mat &
    v()
    {
        return m_v;
    }

    mat &
    F()
    {
        return m_F;
    }

    vec &
    stableMass()
    {
        return m_stableMass;
    }

    mat &
    Fold()
    {
        return m_Fold;
    }

    mat &
    b()
    {
        return m_b;
    }

    mat &
    u()
    {
        return m_u;
    }

    int
    getPdParamId(string paramId)
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
        for(unsigned int col=0;col<m_nParticles; col++)
        {
            for(auto &con:m_PdConnections[col])
            {
                con.second.push_back(value);
            }
        }

        return pos;
    }

    void
    dimensionalScaling(const double E0, const double L0, const double v0,
                       const double t0, const double rho0);
};
//------------------------------------------------------------------------------
}
#endif // PD_PARTICLES_H
