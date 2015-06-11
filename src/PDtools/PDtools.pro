TEMPLATE = lib
TARGET = PDtools
CONFIG += staticlib
CONFIG   -= app_bundle
CONFIG   -= qt

include(../../defaults.pri)

DEFINES *= ARMA_NO_DEBUG
#LIBS *=  -lboost_system -lboost_filesystem

HEADERS += \
    PDtools.h \
    Particles/particles.h \
    Particles/pd_particles.h \
    Grid/grid.h \
    Domain/domain.h \
    Solver/solver.h \
    data/pdshareddata.h \
    Particles/saveparticles.h \
    Particles/loadparticles.h \
    Particles/loadpdparticles.h \
    PdFunctions/pdfunctions.h \
    Solver/timeintegrator.h \
    Solver/TimeIntegrators/velocityverletintegrator.h \
    Solver/solvers.h \
    Force/force.h \
    Force/forces.h \
    Force/PdForces/pd_bondforce.h \
    Modfiers/modifier.h \
    Modfiers/Implementation/velocityboundary.h \
    Modfiers/Implementation/pmbfracture.h \
    SavePdData/savepddata.h \
    SavePdData/Implementations/computedamage.h \
    SavePdData/Implementations/computekineticenergy.h \
    SavePdData/Implementations/computepotentialenergy.h \
    SavePdData/Implementations/computestress.h \
    Force/PdForces/contactforce.h \
    Solver/adr.h \
    Modfiers/Implementation/moveparticles.h \
    Modfiers/modifiers.h \
    Modfiers/Implementation/adrfracture.h \
    SavePdData/Implementations/computemaxstretch.h \
    Modfiers/Implementation/mohrcoulombfracture.h \
    SavePdData/Implementations/computeaveragestretch.h \
    Modfiers/Implementation/adrfractureaverage.h \
    Modfiers/Implementation/adrmohrcoulombfracture.h \
    Force/PdForces/pd_bondforcegaussian.h

SOURCES += \
    Grid/grid.cpp \
    Domain/domain.cpp \
    Solver/solver.cpp \
    data/pdshareddata.cpp \
    Particles/saveparticles.cpp \
    Particles/loadparticles.cpp \
    Particles/loadpdparticles.cpp \
    PdFunctions/pdfunctions.cpp \
    Solver/timeintegrator.cpp \
    Solver/TimeIntegrators/velocityverletintegrator.cpp \
    Force/force.cpp \
    Force/PdForces/pd_bondforce.cpp \
    Modfiers/modifier.cpp \
    Modfiers/Implementation/velocityboundary.cpp \
    Modfiers/Implementation/pmbfracture.cpp \
    SavePdData/savepddata.cpp \
    SavePdData/Implementations/computedamage.cpp \
    SavePdData/Implementations/computekineticenergy.cpp \
    SavePdData/Implementations/computepotentialenergy.cpp \
    SavePdData/Implementations/computestress.cpp \
    Force/PdForces/contactforce.cpp \
    Particles/particles.cpp \
    Solver/adr.cpp \
    Modfiers/Implementation/moveparticles.cpp \
    Modfiers/Implementation/adrfracture.cpp \
    Particles/pd_particles.cpp \
    SavePdData/Implementations/computemaxstretch.cpp \
    Modfiers/Implementation/mohrcoulombfracture.cpp \
    SavePdData/Implementations/computeaveragestretch.cpp \
    Modfiers/Implementation/adrfractureaverage.cpp \
    Modfiers/Implementation/adrmohrcoulombfracture.cpp \
    Force/PdForces/pd_bondforcegaussian.cpp

#headers.path    = $$OUT_PWD
#headers.files   += $$HEADERS
#INSTALLS       += headers

INSTALL_PREFIX = $$OUT_PWD

for(header, HEADERS) {
  path = $${INSTALL_PREFIX}/$${dirname(header)}
  eval(headers_$${path}.files += $$header)
  eval(headers_$${path}.path = $$path)
  eval(INSTALLS *= headers_$${path})
}
