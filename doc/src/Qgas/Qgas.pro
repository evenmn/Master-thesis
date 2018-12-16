TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    hastings_tools.cpp \
    gibbs_tools.cpp \
    test.cpp \
    general_tools.cpp \
    energy.cpp \
    vmc.cpp \
    optimization.cpp \
    basis.cpp

HEADERS += \
    wavefunction.h \
    hastings_tools.h \
    gibbs_tools.h \
    test.h \
    general_tools.h \
    energy.h \
    vmc.h \
    optimization.h \
    basis.h
