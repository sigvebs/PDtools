CONFIG *= OPENMP
#DEFINES *= USE_N3L

contains(CONFIG, OPENMP):message(Building with OPENMP)
contains(DEFINES, USE_N3L):message(Building with N3L)
message(config = $$QMAKESPEC)
#-------------------------------------------------------------------------------
# Optimizations
#-------------------------------------------------------------------------------
CONFIG(debug, debug|release){
    message(Building in debug mode)
    DEFINES *= DEBUG_MODE
} else {
#    DEFINES *= ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
    message(release)
}


CONFIG(OPENMP) {
    message(Building with openMP support)
    DEFINES *= USE_OPENMP

    linux-icc-64 {
        LIBS *= -liomp5 -lpthread
        QMAKE_CXX *= -openmp
    }
    linux-g++ {
        LIBS *= -fopenmp
        QMAKE_CXX *= -fopenmp
    }
}
QMAKE_CXX *= -Wno-unused-result

CONFIG *= c++11
#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
SRC_DIR = $$PWD/src
INCLUDEPATH *= $$SRC_DIR
INCLUDEPATH *= $$SRC_DIR/PDtools

#LIBS *= -lboost_system -lboost_filesystem -lboost_regex
LIBS *= -lboost_regex
LIBS *= -lconfig++
LIBS *= -larmadillo

linux-icc-64 {
#    LIBS *= -lpthread -lmkl_intel_lp64
    LIBS *= -lpthread
}

#-------------------------------------------------------------------------------
# Custom defines
#-------------------------------------------------------------------------------
DEFINES *= PARTICLE_BUFFER=1.3
DEFINES *= PARAMETER_BUFFER=60
DEFINES *= DIM=3
DEFINES *= PD_MAX_CONNECTIONS=120
DEFINES *= PD_MAX_DATA=4
#DEFINES *= ARMA_DONT_USE_WRAPPER
#DEFINES *= ARMA_USE_BLAS
#DEFINES *= ARMA_USE_LAPACK


#INCLUDEPATH *= /usr/include/armadillo_bits
#LIBS *= -larmadillo
#LIBS *= -lconfig++
#LIBS *= -lgfortran -llapack -lopenblas
#LIBS *= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#LIBS += -L$$TOP_OUT_PWD/src/PDtools -lPDtools
#QMAKE_CXXFLAGS += -vec-report=5
#QMAKE_LFLAGS += -vec-report=5


