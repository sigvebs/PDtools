#ifndef SAVESTATE_H
#define SAVESTATE_H



#include "PDtools/Modfiers/modifier.h"

namespace PDtools
{
//------------------------------------------------------------------------------
class SaveState : public Modifier
{
public:
    SaveState(int frequency);

    virtual void
    initialize();

    virtual void
    evaluateStepOne();
protected:
    int m_frequency;
    int m_counter;
};
//------------------------------------------------------------------------------
}
#endif // SAVESTATE_H
