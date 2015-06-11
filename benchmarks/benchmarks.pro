TEMPLATE  = app
TARGET 	  = tests
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG   += thread

include(../defaults.pri)
DEFINES *= ARMA_NO_DEBUG

LIBS += $$TOP_OUT_PWD/src/PDtools/libPDtools.a
LIBS += -lgtest

SOURCES += \
    main.cpp \
    armadillo/testarmadillo.cpp

HEADERS += \
    test_resources.h
