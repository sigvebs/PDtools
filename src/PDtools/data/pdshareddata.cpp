#include "pdshareddata.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PDsharedData::PDsharedData():
    m_domain(NULL),
    m_pdGrid(NULL),
    m_particles(NULL)
{

}
//------------------------------------------------------------------------------
PDsharedData::~PDsharedData()
{
    if(m_domain != NULL)
    {
        delete m_domain;
    }
    if(m_pdGrid != NULL)
    {
        delete m_pdGrid;
    }
    if(m_particles != NULL)
    {
        delete m_particles;
    }
}
//------------------------------------------------------------------------------
void PDsharedData::registerNew(const string &key, double value)
{
    assert(m_data.count(key) == 0);
    m_data[key] = value;
}
//------------------------------------------------------------------------------
void PDsharedData::updateData(const string &key, double value)
{
    assert(m_data.count(key) == 1);
    m_data[key] = value;
}
//------------------------------------------------------------------------------
}
