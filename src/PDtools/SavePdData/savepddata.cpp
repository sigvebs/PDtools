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
            m_saveparam_scale.push_back(std::pair<std::string, double>("damage", 1.));
        }
        if(param == "max_stretch")
        {
            m_computeProperties.push_back(new ComputeMaxStretch(*m_particles));
            m_saveparam_scale.push_back(std::pair<std::string, double>("max_stretch", 1.));
        }
        if(param == "radius")
        {
            m_saveparam_scale.push_back(std::pair<std::string, double>("radius", m_L0));
        }
        if(param == "average_stretch")
        {
            m_computeProperties.push_back(new ComputeAverageStretch(*m_particles));
            m_saveparam_scale.push_back(std::pair<std::string, double>("average_stretch", 1.));
        }
        else if(param == "kinetic_energy")
        {
            m_computeProperties.push_back(new ComputeKineticEnergy(*m_particles));
            m_saveparam_scale.push_back(std::pair<std::string, double>("kinetic_energy", pow(m_v0, 2)*pow(m_L0, 3)*m_rho0));
        }
        else if(param == "potential_energy")
        {
            m_computeProperties.push_back(new ComputePotentialEnergy(*m_particles,
                                                                     *m_oneBodyForces));
            m_saveparam_scale.push_back(std::pair<std::string, double>("potential_energy", 1.));
        }
        else if(param == "stress")
        {
            m_computeProperties.push_back(new ComputeStress(*m_particles,
                                                            *m_oneBodyForces));
            computeStress = true;
        }
        else if(param == "id")
        {
            m_saveparam_scale.push_back(std::pair<std::string, double>("id", 1.));
        }
        else if(param == "x" || param == "y" || param == "z")
        {
            m_saveparam_scale.push_back(std::pair<std::string, double>(param, m_L0));
        }
        else if(   param == "s_xx" || param == "s_yy" || param == "s_zz"
                || param == "s_xy" || param == "s_xz" || param == "s_yz")
        {
            m_saveparam_scale.push_back(std::pair<std::string, double>(param, m_E0));
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
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_xx", m_E0));
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_yy", m_E0));
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_zz", m_E0));
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_xy", m_E0));
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_xz", m_E0));
        m_saveparam_scale.push_back(std::pair<std::string, double>("s_yz", m_E0));
    }
#ifdef USE_OPENMP
# pragma omp parallel for
#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        pair<int, int> id(i, i);
        for(ComputeProperty *computeProperty:m_computeProperties)
        {
            computeProperty->init(id);
        }
    }

    saveParticles = new SaveParticles("lmp", m_saveparam_scale, false);
}
//------------------------------------------------------------------------------
void SavePdData::evaluate(double t, int i)
{
    (void) t;

//#ifdef USE_OPENMP
//# pragma omp parallel for
//#endif
    for(unsigned int i=0; i<m_particles->nParticles(); i++)
    {
        pair<int, int> id(i, i);
        for(ComputeProperty *computeProperty:m_computeProperties)
        {
            computeProperty->update(id);
        }
    }
}
//------------------------------------------------------------------------------
void SavePdData::saveData(double t, int i)
{
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
void SavePdData::setScaling(const double E0, const double L0, const double v0, const double t0, const double rho0)
{
    m_E0 = E0;
    m_L0 = L0;
    m_v0 = v0;
    m_t0 = t0;
    m_rho0 = rho0;
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
