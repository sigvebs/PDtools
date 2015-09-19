#include "savepddata.h"

#include "PDtools/SavePdData/Implementations/computedamage.h"
#include "PDtools/SavePdData/Implementations/computekineticenergy.h"
#include "PDtools/SavePdData/Implementations/computepotentialenergy.h"
#include "PDtools/SavePdData/Implementations/computestress.h"
#include "PDtools/SavePdData/Implementations/computemaxstretch.h"
#include "PDtools/SavePdData/Implementations/computeaveragestretch.h"

namespace PDtools
{
//------------------------------------------------------------------------------
SavePdData::SavePdData()
{

}
//------------------------------------------------------------------------------
SavePdData::SavePdData(std::vector<string> saveParameters):
    SavePdData()
{
    m_saveParameters = saveParameters;
}
//------------------------------------------------------------------------------
SavePdData::~SavePdData()
{
    // delete m_computeProperties
    delete saveParticles;
}
//------------------------------------------------------------------------------
void SavePdData::initialize()
{
    if(m_particles == nullptr)
    {
        cerr << "Particles not set for savePdData" << endl;
        throw;
    }

//    boost::filesystem::path dir(m_savePath);
//    if(boost::filesystem::create_directories(dir)) {
//        std::cout << "Directory created: " << m_savePath << "\n";
//    }

    bool computeStress = false;

    for(string param:m_saveParameters)
    {
        if(param == "damage")
        {
            m_computeProperties.push_back(new ComputeDamage(*m_particles));
        }
        if(param == "max_stretch")
        {
            m_computeProperties.push_back(new ComputeMaxStretch(*m_particles));
        }
        if(param == "f_x")
        {
            m_computeProperties.push_back(new ComputeMaxStretch(*m_particles));
        }
        if(param == "average_stretch")
        {
            m_computeProperties.push_back(new ComputeAverageStretch(*m_particles));
        }
        else if(param == "kinetic_energy")
        {
            m_computeProperties.push_back(new ComputeKineticEnergy(*m_particles));
        }
        else if(param == "potential_energy")
        {
            m_computeProperties.push_back(new ComputePotentialEnergy(*m_particles,
                                                                     *m_oneBodyForces));
        }
        else if(param == "stress")
        {
            m_computeProperties.push_back(new ComputeStress(*m_particles,
                                                            *m_oneBodyForces));
            computeStress = true;
        }
        else if(param == "id"
                || param == "x" || param == "y" || param == "z"
//                || param == "v_x" || param == "v_y" || param == "v_z"
                || param == "s_xx" || param == "s_yy" || param == "s_zz"
                || param == "s_xy" || param == "s_xz" || param == "s_yz"
                )
        {
            // Ignore id, velocities, coordinates and stress
        }
        else
        {
            if(!m_particles->hasParameter(param))
            {
                cerr << "ERROR: In savePdData. Particles does not containt the "
                     << "parameter '" << param << "' and param is not a "
                     << "calculabel type.";
                throw;
            }
        }
    }

    // Adding the stress parameters
    if(computeStress)
    {
        m_saveParameters.push_back("s_xx");
        m_saveParameters.push_back("s_yy");
        m_saveParameters.push_back("s_zz");
        m_saveParameters.push_back("s_xy");
        m_saveParameters.push_back("s_xz");
        m_saveParameters.push_back("s_yz");
    }
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        pair<int, int> id(i, i);
        for(ComputeProperty *computeProperty:m_computeProperties)
        {
            computeProperty->init(id);
        }
    }

    saveParticles = new SaveParticles("lmp", m_saveParameters, false);
}
//------------------------------------------------------------------------------
void SavePdData::evaluate(double t, int i)
{
    (void) t;

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(int i=0; i<m_particles->nParticles(); i++)
    {
        pair<int, int> id(i, i);
        for(ComputeProperty *computeProperty:m_computeProperties)
        {
            computeProperty->update(id);
        }
    }

    string fName = m_savePath + "/" + to_string(i) + ".lmp";
    saveParticles->setTimestep(i);
    saveParticles->writeToFile(*m_particles, fName);
}
//------------------------------------------------------------------------------
void SavePdData::setSavePath(string savePath)
{
    m_savePath = savePath;
}
//------------------------------------------------------------------------------
void SavePdData::setParticles(PD_Particles *particles)
{
    m_particles = particles;
}
//------------------------------------------------------------------------------
void SavePdData::addParameter(string param)
{
    m_saveParameters.push_back(param);
}
//------------------------------------------------------------------------------
void SavePdData::setForces(std::vector<Force *> &oneBodyForces)
{
    m_oneBodyForces = &oneBodyForces;
}
//------------------------------------------------------------------------------
// Compute Property class
//------------------------------------------------------------------------------
ComputeProperty::ComputeProperty(PD_Particles &particles):
    m_particles(particles)
{

}
//------------------------------------------------------------------------------
void ComputeProperty::init(const pair<int, int> &pIdcol)
{
    (void) pIdcol;
}
//------------------------------------------------------------------------------
}
