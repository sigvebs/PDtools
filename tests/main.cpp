#include <PDtools.h>
#include <vector>
#include <armadillo>
#include <stdio.h>
#include <gtest/gtest.h>
#include "test_resources.h"

using namespace std;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
