TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    InitialStates/randomnormal.cpp \
    WaveFunctions/hydrogenorbital.cpp \
    Hamiltonians/atomicnucleus.cpp \
    WaveFunctions/padejastrow.cpp \
    Math/random2.cpp \
    InitialWeights/initialweights.cpp \
    InitialWeights/ones.cpp \
    InitialWeights/randomize.cpp \
    WaveFunctions/mlgaussian.cpp \
    WaveFunctions/padejastrowcartesian.cpp \
    WaveFunctions/cartesiangaussian.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    InitialStates/randomnormal.h \
    WaveFunctions/hydrogenorbital.h \
    Hamiltonians/atomicnucleus.h \
    WaveFunctions/padejastrow.h \
    Math/random2.h \
    InitialWeights/initialweights.h \
    InitialWeights/ones.h \
    InitialWeights/randomize.h \
    WaveFunctions/mlgaussian.h \
    WaveFunctions/padejastrowcartesian.h \
    WaveFunctions/cartesiangaussian.h

