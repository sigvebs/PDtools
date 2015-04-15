#include "pdshareddata.h"

namespace PDtools
{
//------------------------------------------------------------------------------
PDsharedData::PDsharedData():
    domain(NULL),
    pdGrid(NULL),
    particles(NULL)
{

}
//------------------------------------------------------------------------------
PDsharedData::~PDsharedData()
{
    if(domain != NULL)
    {
        delete domain;
    }
    if(pdGrid != NULL)
    {
        delete pdGrid;
    }
    if(particles != NULL)
    {
        delete particles;
    }
}
//------------------------------------------------------------------------------
void PDsharedData::registerNew(const string &key, double value)
{
    assert(m_data.count(key) == 0);    m_data[key] = value;
}
//------------------------------------------------------------------------------
void PDsharedData::updateData(const string &key, double value)
{
    assert(m_data.count(key) == 1);
    m_data[key] = value;
}
//------------------------------------------------------------------------------
}
