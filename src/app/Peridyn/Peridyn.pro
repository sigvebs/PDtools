TEMPLATE  = app
TARGET 	  = peridyn
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

INCLUDEPATH += $$SRC_DIR/PDtools
LIBS *= $$TOP_OUT_PWD/src/PDtools/libPDtools.a
#LIBS *= -lboost_regex
#LIBS *= -lconfig++
include(../../../defaults.pri)

# For dynamic linking:
#LIBS +=  -L$$TOP_OUT_PWD/src/PDtools -lPDtools
#QMAKE_LFLAGS += -Wl,--rpath=$$TOP_OUT_PWD/src/PDtools

#message(from peridyn.pro $$LIBS)

SOURCES += \
    main.cpp \
    pdsolver.cpp

HEADERS += \
    pdsolver.h
