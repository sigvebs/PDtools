CONFIG += OPENMP
#DEFINES *= USE_N3L

#-------------------------------------------------------------------------------
# Optimizations
#-------------------------------------------------------------------------------

CONFIG(debug, debug|release){
    message(Building in debug mode)
    DEFINES *= DEBUG_MODE
} else {
    DEFINES *= ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
    message(release)
}

CONFIG(OPENMP) {
    message(Building with openMP support)
    DEFINES *= USE_OPENMP
    LIBS += -fopenmp
    QMAKE_CXX += -fopenmp
}

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
SRC_DIR = $$PWD/src
INCLUDEPATH += $$SRC_DIR
INCLUDEPATH += $$SRC_DIR/PDtools

# For python includes - should be general
#INCLUDEPATH += -I /usr/include/python2.7/
#LIBS *= -lboost_program_options -lboost_python

LIBS *=  -lconfig++ -larmadillo -llapack -lblas
LIBS *=  -lboost_system -lboost_filesystem

COMMON_CXXFLAGS += -std=c++11
COMMON_CXXFLAGS += -Wno-sign-compare

QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

#LIBS += -L$$TOP_OUT_PWD/src/PDtools -lPDtools

#QMAKE_CXXFLAGS += -vec-report=5
#QMAKE_LFLAGS += -vec-report=5

#-------------------------------------------------------------------------------
# Custom defines
#-------------------------------------------------------------------------------
DEFINES += PARTICLE_BUFFER=1.3 \
    PARAMETER_BUFFER=60 \
    DIM=3 \
    PD_MAX_CONNECTIONS=120 \
    PD_MAX_DATA=4
