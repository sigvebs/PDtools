#include "calculatepdangles.h"

#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculatePdAngles::CalculatePdAngles():
    CalculateProperty("PdAngle")
{

}
//------------------------------------------------------------------------------
CalculatePdAngles::~CalculatePdAngles()
{

}
//------------------------------------------------------------------------------
void CalculatePdAngles::initialize()
{
    m_indexConnected = m_particles->getPdParamId("connected");

    if(m_dim==2)
    {
        m_indexTheta = m_particles->registerPdParameter("theta");
        m_indexTheta0 = m_particles->registerPdParameter("theta0");
    }
    else
    {
        cerr << "Perperties - calculatePdAngles: dim " << m_dim
             << " is not supported"
             << endl;
    }

    // Calculating the initial angle
    const int nParticles = m_particles->nParticles();
    std::unordered_map<int, int> & idToCol = m_particles->idToCol();
    const arma::ivec & colToId= m_particles->colToId();

    const mat & R = m_particles->r0();

    double dr_ij[m_dim];
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId[i];
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
        const int nConnections = PDconnections.size();

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            auto &con = PDconnections[l_j];
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = idToCol.at(id_j);

            for(int d=0;d<m_dim;d++)
            {
                dr_ij[d] = R(j, d) - R(i, d);
            }
            const double theta = atan2(dr_ij[1], dr_ij[0]);
            con.second[m_indexTheta0] = theta;
        }
    }
}
//------------------------------------------------------------------------------
void CalculatePdAngles::update()
{
    const int nParticles = m_particles->nParticles();
    std::unordered_map<int, int> & idToCol = m_particles->idToCol();
    const arma::ivec & colToId= m_particles->colToId();

    const mat & R = m_particles->r();

    double dr_ij[m_dim];
    for(int i=0; i<nParticles; i++)
    {
        const int id_i = colToId[i];
        vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
        const int nConnections = PDconnections.size();

        for(int l_j=0; l_j<nConnections; l_j++)
        {
            auto &con = PDconnections[l_j];
            if(con.second[m_indexConnected] <= 0.5)
                continue;

            const int id_j = con.first;
            const int j = idToCol.at(id_j);

            for(int d=0;d<m_dim;d++)
            {
                dr_ij[d] = R(j, d) - R(i, d);
            }
            const double theta = atan2(dr_ij[1], dr_ij[0]);
            con.second[m_indexTheta] = theta - con.second[m_indexTheta0];
        }
    }
}
//------------------------------------------------------------------------------
}

