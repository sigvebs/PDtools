#include <gtest/gtest.h>
#include <PDtools.h>

using namespace PDtools;

extern std::vector<string> geometries;


class GRID_FIXTURE : public ::testing::Test {
    protected:
    Domain *domain;
    std::vector<std::pair<double,double>> boundaries;

    GRID_FIXTURE()
    {
        std::pair<double,double> x_limits(-5., 2.1);
        std::pair<double,double> y_limits(0, 3.6);
        std::pair<double,double> z_limits(-1., 1.2);
        boundaries = {x_limits, y_limits, z_limits};
        domain = new Domain(3, boundaries);
    }

    ~GRID_FIXTURE()
    {
        delete domain;
    }
};


TEST_F(GRID_FIXTURE, PD_DOMAIN)
{
    ASSERT_EQ(domain->dim, 3)
            << "The domain dimension is wrong.";
    const vector<std::pair<double,double>> boundaries = domain->boundaries();
    ASSERT_EQ(boundaries[0].first, -5)
            << "The domain limits are wrong.";
    ASSERT_EQ(boundaries[0].second, 2.1)
            << "The domain limits are wrong.";
    ASSERT_EQ(boundaries[1].first, 0)
            << "The domain limits are wrong.";
    ASSERT_EQ(boundaries[1].second, 3.6)
            << "The domain limits are wrong.";
    ASSERT_EQ(boundaries[2].first, -1)
            << "The domain limits are wrong.";
    ASSERT_EQ(boundaries[2].second, 1.2)
            << "The domain limits are wrong.";
}


TEST_F(GRID_FIXTURE, PD_GRID)
{
//    arma::ivec test_innerGripoints_1 = {0, 4, 1, 5, 2, 6, 3, 7};
//    arma::ivec test_innerGripoints_2 = {4, 8, 5, 9, 6, 10, 7, 11};
//    arma::ivec test_ghostGripoints_2 = {0, 1, 2, 3, 12, 13, 14, 15};

//    double gridSpacing = 1.5;
//    Grid grid(boundaries, gridSpacing);
//    grid.createGrid();

//    unordered_map<int, GridPoint> & gridpoints_1 = grid.gridpoints();

//    ASSERT_EQ(gridpoints_1.size(), 8);

//    for(int testPoint:test_innerGripoints_1)
//        ASSERT_EQ(gridpoints_1.count(testPoint), 1);

//    arma::ivec3 pb = {0,1,0};
//    Grid grid_pb(boundaries, gridSpacing, pb);
//    grid_pb.update();
//    grid_pb.createGrid();

//    unordered_map<int, GridPoint> & gridpoints_2 = grid_pb.gridpoints();

//    ASSERT_EQ(gridpoints_2.size(), 16);
//    for(int testPoint:test_innerGripoints_2)
//    {
//        ASSERT_EQ(gridpoints_2.count(testPoint), 1);
//    }

//    for(int testPoint:test_ghostGripoints_2)
//    {
//        ASSERT_EQ(gridpoints_2.count(testPoint), 1);
//        ASSERT_EQ(gridpoints_2[testPoint].isGhost(), 1);
//    }
}

TEST_F(GRID_FIXTURE, PD_GRID_DOMAIN)
{
    arma::ivec test_innerGripoints_1 = {0, 4, 1, 5, 2, 6, 3, 7};
    arma::ivec test_innerGripoints_2 = {4, 8, 5, 9, 6, 10, 7, 11};
    arma::ivec test_ghostGripoints_2 = {0, 1, 2, 3, 12, 13, 14, 15};

    double gridSpacing = 1.5;
    Grid grid(*domain, gridSpacing);

    grid.createGrid();
//    unordered_map<int, GridPoint> & gridpoints_1 = grid.gridpoints();

//    ASSERT_EQ(gridpoints_1.size(), 8);

//    for(int testPoint:test_innerGripoints_1)
//        ASSERT_EQ(gridpoints_1.count(testPoint), 1);

    arma::ivec3 pb = {0,1,0};
    domain->periodicBoundaries(pb);
    Grid grid_pb(*domain, gridSpacing);
    grid_pb.createGrid();

    unordered_map<int, GridPoint*> & gridpoints_2 = grid_pb.gridpoints();

//    ASSERT_EQ(gridpoints_2.size(), 16);
//    for(int testPoint:test_innerGripoints_2)
//    {
//        ASSERT_EQ(gridpoints_2.count(testPoint), 1);
//    }

//    for(int testPoint:test_ghostGripoints_2)
//    {
//        ASSERT_EQ(gridpoints_2.count(testPoint), 1);
//        ASSERT_EQ(gridpoints_2[testPoint].isGhost(), 1);
//    }
}


TEST_F(GRID_FIXTURE, PD_GRID_NEIGHBOURS)
{
    unordered_map<int, arma::ivec> test_neigbours;

    test_neigbours[4] = arma::ivec({0, 8, 1, 5, 9});
    test_neigbours[8] = arma::ivec({4, 12, 5, 9, 13});
    test_neigbours[5] = arma::ivec({0, 4, 8, 1, 9, 2, 6, 10});
    test_neigbours[9] = arma::ivec({4, 8, 12, 5, 13, 6, 10, 14});
    test_neigbours[6] = arma::ivec({1, 5, 9, 2, 10, 3, 7, 11});
    test_neigbours[10] = arma::ivec({5, 9, 13, 6, 14, 7, 11, 15});
    test_neigbours[7] = arma::ivec({2, 6, 10, 3, 11});
    test_neigbours[11] = arma::ivec({6, 10, 14, 7, 15});

    arma::ivec3 pb = {0,1,0};
    double gridSpacing = 1.5;
    domain->periodicBoundaries(pb);
    Grid grid(*domain, gridSpacing);
    grid.initialize();

    unordered_map<int, GridPoint*> & gridpoints = grid.gridpoints();

    for(pair<int, arma::ivec> test_idNeighbours:test_neigbours)
    {
//        int id = test_idNeighbours.first;
//        ASSERT_EQ(gridpoints.count(id), 1);
//        vector<pair<int, GridPoint&>> neighbours = gridpoints[id].neighbours();

//        for(int i=0; i<neighbours.size(); i++)
//        {
//            ASSERT_EQ(neighbours[i].first, test_neigbours[id][i]);
//            ++i;
//        }
    }
}
