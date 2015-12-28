#include <iostream>
#include <armadillo>
#include <memory>

#include "PDtools.h"

using namespace std;

int main(int argc, char** argv)
{
    (void) argv;
    if(argc != 2)
    {
        cerr << "usage: peridyn [path to config file]" << endl;
        return 1;
    }

     return 0;
}
