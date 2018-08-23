TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    gradient_descent.cpp \
    hastings_tools.cpp \
    gibbs_tools.cpp \
    test.cpp \
    basis.cpp

HEADERS += \
    wavefunction.h \
    gradient_descent.h \
    hastings_tools.h \
    gibbs_tools.h \
    test.h \
    basis.h
