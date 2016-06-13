cache()

TEMPLATE = subdirs
SUBDIRS = src
#tests benchmarks
#tests.depends = src
CONFIG += ordered

src.subdirs = pd_lib
test.depends = pd_lib
#benchmarks.depends = pd_lib

OTHER_FILES += \
    defaults.pri

