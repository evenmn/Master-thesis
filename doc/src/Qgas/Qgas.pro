TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    hastings_tools.cpp \
    gibbs_tools.cpp \
    test.cpp \
    basis.cpp \
    general_tools.cpp \
    VMC.cpp \
    gradient_descent.cpp

HEADERS += \
    wavefunction.h \
    hastings_tools.h \
    gibbs_tools.h \
    test.h \
    basis.h \
    general_tools.h \
    VMC.h \
    gradient_descent.h
