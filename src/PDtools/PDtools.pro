TEMPLATE = lib
TARGET = PDtools
CONFIG += staticlib
CONFIG   -= app_bundle
CONFIG   -= qt

include(../../defaults.pri)

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
    Force/PdForces/pd_bondforcegaussian.h \
    Force/PdForces/pd_pmb.h \
    Force/PdForces/pd_osp.h \
    Modfiers/Implementation/boundaryforce.h \
    Modfiers/Implementation/bondenergyfracture.h \
    Modfiers/Implementation/simplefracture.h \
    Force/PdForces/pd_lps.h \
    Solver/ADRsolvers/dynamicadr.h \
    Solver/staticsolver.h \
    Solver/TimeIntegrators/eulercromerintegrator.h \
    PdFunctions/pdfunctionsmpi.h \
    SavePdData/Implementations/computegridid.h \
    Force/PdForces/viscousdamper.h \
    CalculateProperties/calculateproperty.h \
    CalculateProperties/Implementation/calculatepdangles.h \
    Modfiers/Implementation/micropolarfracture.h \
    CalculateProperties/calculateproperties.h \
    Force/DemForces/demforce.h \
    CalculateProperties/Implementation/calculatestress.h \
    Modfiers/Implementation/mohrcoulombbondfracture.h \
    Modfiers/Implementation/rigidwall.h \
    Modfiers/Implementation/adrmohrcoulombbondfracture.h \
    Modfiers/Implementation/mohrcoulombmax.h \
    Modfiers/Implementation/mohrcoulombmaxfracture.h \
    Modfiers/Implementation/mohrcoulombweightedaverage.h \
    Force/PdForces/pd_dampenedbondforce.h

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
    Force/PdForces/pd_bondforcegaussian.cpp \
    Force/PdForces/pd_pmb.cpp \
    Force/PdForces/pd_osp.cpp \
    Modfiers/Implementation/boundaryforce.cpp \
    Modfiers/Implementation/bondenergyfracture.cpp \
    Modfiers/Implementation/simplefracture.cpp \
    Force/PdForces/pd_lps.cpp \
    Solver/ADRsolvers/dynamicadr.cpp \
    Solver/staticsolver.cpp \
    Solver/TimeIntegrators/eulercromerintegrator.cpp \
    PdFunctions/pdfunctionsmpi.cpp \
    SavePdData/Implementations/computegridid.cpp \
    Force/PdForces/viscousdamper.cpp \
    CalculateProperties/calculateproperty.cpp \
    CalculateProperties/Implementation/calculatepdangles.cpp \
    Modfiers/Implementation/micropolarfracture.cpp \
    Force/DemForces/demforce.cpp \
    CalculateProperties/Implementation/calculatestress.cpp \
    Modfiers/Implementation/mohrcoulombbondfracture.cpp \
    Modfiers/Implementation/rigidwall.cpp \
    Modfiers/Implementation/adrmohrcoulombbondfracture.cpp \
    Modfiers/Implementation/mohrcoulombmax.cpp \
    Modfiers/Implementation/mohrcoulombmaxfracture.cpp \
    Modfiers/Implementation/mohrcoulombweightedaverage.cpp \
    Force/PdForces/pd_dampenedbondforce.cpp

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
