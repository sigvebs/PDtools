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
//    ::testing::GTEST_FLAG(filter) = "PD_SOLVER_FIXTURE*";
    ::testing::GTEST_FLAG(filter) = "PD_STRETCH_FIXTURE*";
    return RUN_ALL_TESTS();
}
