TEMPLATE  = app
TARGET 	  = tests
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG   += thread

include(../defaults.pri)
#DEFINES *= ARMA_NO_DEBUG

LIBS *= $$TOP_OUT_PWD/src/PDtools/libPDtools.a
LIBS *= -lgtest

SOURCES += \
    main.cpp \
#    PDtools/particles/test_particles.cpp \
#    PDtools/PD_particles/test_pd_particles.cpp \
#    PDtools/grid/test_grid.cpp \
#    PDtools/configuration/test_configuration/test_configuration.cpp \
#    PDtools/test_solver/test_solver.cpp \
#    PDtools/stretch/test_stretch.cpp \
#    PDtools/LinearSolver/linearsolver.cpp \
#    PDtools/MPI/test_mpi.cpp

#HEADERS += \
#    test_resources.h

#-------------------------------------------------------------------------------
# Creates extra test directories
#-------------------------------------------------------------------------------
MKDIR = /bin/mkdir
DIRECTORIES = testGeometries
DEFINES += TEST_SAVE_PATH=\\\"testGeometries\\\"

for(DIRECTORY, DIRECTORIES) {
     mkcommands += $$OUT_PWD/$$DIRECTORY
}

createDirs.commands = $(MKDIR) $$mkcommands
first.depends += createDirs
QMAKE_EXTRA_TARGETS += first createDirs

message($$DEFINES)
#-------------------------------------------------------------------------------
# Compiling the google test framework on Ubuntu
#-------------------------------------------------------------------------------
# sudo apt-get install libgtest-dev
#cd /usr/src/gtest
#sudo cmake .
#sudo make
#sudo mv libg* /usr/lib/
