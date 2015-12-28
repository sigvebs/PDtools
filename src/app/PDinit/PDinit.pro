TEMPLATE  = app
TARGET 	  = PDinit
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

include(../../../defaults.pri)

INCLUDEPATH += $$SRC_DIR/PDtools
LIBS += -lconfig++
LIBS += $$TOP_OUT_PWD/src/PDtools/libPDtools.a

SOURCES += \
    main.cpp

