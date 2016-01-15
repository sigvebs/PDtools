#CONFIG *= OPENMP
CONFIG *= MPI
#DEFINES *= USE_N3L
#DEFINES *= USE_BOOST_MPI

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
    message(Building in release mode)
    DEFINES *= ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
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

CONFIG(abel) {
    INCLUDEPATH += /cluster/home/sigve/usr/local/include
    INCLUDEPATH += /usit/abel/u1/sigve/usr/local/include
    INCLUDEPATH += /usit/abel/u1/sigve/usr/boost/include

    INCLUDEPATH += /usit/abel/u1/sigve/usr/include
    INCLUDEPATH += /usit/abel/u1/sigve/usr/local/include
    INCLUDEPATH += /usit/abel/u1/sigve/usr/boost/include
    LIBS *= -L/usit/abel/u1/sigve/usr/local/lib -lconfig++
    LIBS *= -L/usit/abel/u1/sigve/usr/lib64 -larmadillo

    QMAKE_CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE += -DMPICH_IGNORE_CXX_SEEK
}

CONFIG(MPI) {
    message(Building with MPI support)
    DEFINES *= USE_MPI

    linux-icc-64 {
        DEFINES -= USE_MPI
        CONFIG -= MPI
#        QMAKE_CXX = mpiicpc
#        QMAKE_CC = mpiicc
#        QMAKE_CXX_RELEASE = $$QMAKE_CXX
#        QMAKE_CXX_DEBUG = $$QMAKE_CXX
#        QMAKE_LINK = $$QMAKE_CXX
#        LIBS *= -parallel

#        CONFIG(debug, debug|release){
#            QMAKE_CXXFLAGS += -g
#        } else {
#            QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

#            QMAKE_CXXFLAGS_RELEASE -= -O1
#            QMAKE_CXXFLAGS_RELEASE -= -O2
#            QMAKE_CXXFLAGS_RELEASE *= -O3
#            QMAKE_CXXFLAGS_RELEASE *= -fast
#            #QMAKE_CXXFLAGS_RELEASE *= -ffast-math
#            #QMAKE_CXXFLAGS_RELEASE *= -funroll-loops
#        }
    }

    linux-g++ {
        QMAKE_CXX = mpicxx
        QMAKE_CC = mpicc
        QMAKE_CXX_RELEASE = $$QMAKE_CXX
        QMAKE_CXX_DEBUG = $$QMAKE_CXX
        QMAKE_LINK = $$QMAKE_CXX

        QMAKE_CFLAGS = $$system(mpicxx --showme:compile)
        QMAKE_LFLAGS = $$system(mpicxx --showme:link)
        QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile)

        CONFIG(debug, debug|release){
            QMAKE_CXXFLAGS += -g
        } else {
            QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

            QMAKE_CXXFLAGS_RELEASE -= -O1
            QMAKE_CXXFLAGS_RELEASE -= -O2
            QMAKE_CXXFLAGS_RELEASE *= -O3
            QMAKE_CXXFLAGS_RELEASE *= -fno-strict-aliasing
            QMAKE_CXXFLAGS_RELEASE *= -ffast-math
            QMAKE_CXXFLAGS_RELEASE *= -funroll-loops
        }
    }
}

CONFIG(MPI) {
    LIBS += -lmpi
    LIBS *= -lboost_system -lboost_filesystem -lboost_serialization

}
linux-g++ {
    QMAKE_CXX *= -Wno-unused-result
}

CONFIG *= c++11
#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
SRC_DIR = $$PWD/src
INCLUDEPATH *= $$SRC_DIR
INCLUDEPATH *= $$SRC_DIR/PDtools

LIBS *= -lboost_regex
LIBS *= -lconfig++
LIBS *= -larmadillo

linux-icc-64 {
#    LIBS *= -lmkl_intel_lp64
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
