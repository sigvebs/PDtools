#include "calculatedamage.h"

#include "Particles/pd_particles.h"

namespace PDtools
{
//------------------------------------------------------------------------------
CalculateDamage::CalculateDamage(double delta):
    CalculateProperty("damage"),
    m_delta(delta)
{

}
//------------------------------------------------------------------------------
void CalculateDamage::initialize()
{
    m_iDamage = m_particles->registerParameter("damage");
    m_iInitialWeight = m_particles->registerParameter("initialWeight");
    m_indexConnected = m_particles->getPdParamId("connected");
    m_iVolumeScaling = m_particles->getPdParamId("volumeScaling");
    m_iVolume = m_particles->getParamId("volume");
    m_iDr0 = m_particles->getPdParamId("dr0");

    const ivec &colToId = m_particles->colToId();
    const unordered_map<int, int> & idToCol = m_particles->idToCol();
    const int nParticles = m_particles->nParticles();
    mat & data = m_particles->data();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nParticles; i++) {
        const int id_i = colToId(i);
        data(i, m_iInitialWeight) = m_particles->pdConnections(id_i).size();

        const vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);

        double initialWeight = 0;
        for (const auto& con:PDconnections) {

            const int id_j = con.first;
            const int j = idToCol.at(id_j);

            const double dr0 = con.second[m_iDr0];
            const double vol_j = data(j, m_iVolume);
            const double volumeScaling = con.second[m_iVolumeScaling];
            const double w = weightFunction(dr0);
//            initialWeight += w*vol_j*volumeScaling;
            initialWeight += vol_j*volumeScaling;
        }
        if(initialWeight == 0)
            initialWeight = 1.0;

        data(i, m_iInitialWeight) = initialWeight;
    }
}
//------------------------------------------------------------------------------
void CalculateDamage::update()
{
    const ivec &colToId = m_particles->colToId();
    const int nParticles = m_particles->nParticles();
    mat & data = m_particles->data();
    const unordered_map<int, int> & idToCol = m_particles->idToCol();

    // Updating single particle states
    for(int i=0; i<nParticles; i++) {
        const int id_i = colToId(i);
        const vector<pair<int, vector<double>>> & PDconnections = m_particles->pdConnections(id_i);
        const int jnum = PDconnections.size();
        const double maxConnections = jnum;

        if(maxConnections <= 0) {
            data(i, m_iDamage) = 1.;
            return;
        }

        double over = 0;
        double under =  data(i, m_iInitialWeight);

        for (const auto& con:PDconnections) {
            if(con.second[m_indexConnected] > 0.5) {
                const int id_j = con.first;
                const int j = idToCol.at(id_j);
                const double dr0 = con.second[m_iDr0];
                const double vol_j = data(j, m_iVolume);
                const double volumeScaling = con.second[m_iVolumeScaling];
                const double w = weightFunction(dr0);

//                over += w*vol_j*volumeScaling;
                over += vol_j*volumeScaling;
            }
        }

        data(i, m_iDamage) = 1. - over/under;
//        data(i, m_indexDamage) = 1. - tot/maxConnections;
        //    data(i, m_indexDamage) = over;
    }
}
//------------------------------------------------------------------------------
}
