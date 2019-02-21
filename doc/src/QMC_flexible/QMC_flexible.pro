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
    InitialStates/randomnormal.cpp \
    WaveFunctions/hydrogenorbital.cpp \
    Hamiltonians/atomicnucleus.cpp \
    Math/random2.cpp \
    InitialWeights/initialweights.cpp \
    InitialWeights/randomize.cpp \
    WaveFunctions/mlgaussian.cpp \
    Metropolis/metropolis.cpp \
    Metropolis/bruteforce.cpp \
    Metropolis/importancesampling.cpp \
    WaveFunctions/nqsjastrow.cpp \
    InitialWeights/constant.cpp \
    Optimization/optimization.cpp \
    Optimization/gradientdescent.cpp \
    WaveFunctions/gaussian.cpp \
    WaveFunctions/padejastrow.cpp \
    WaveFunctions/partlyrestricted.cpp \
    WaveFunctions/slaterdeterminant.cpp \
    Optimization/barzilaiborwein.cpp \
    Optimization/gradientdescentmomentum.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    InitialStates/randomnormal.h \
    WaveFunctions/hydrogenorbital.h \
    Hamiltonians/atomicnucleus.h \
    Math/random2.h \
    InitialWeights/initialweights.h \
    InitialWeights/randomize.h \
    WaveFunctions/mlgaussian.h \
    Metropolis/metropolis.h \
    Metropolis/bruteforce.h \
    Metropolis/importancesampling.h \
    WaveFunctions/nqsjastrow.h \
    InitialWeights/constant.h \
    Optimization/optimization.h \
    Optimization/gradientdescent.h \
    WaveFunctions/gaussian.h \
    WaveFunctions/padejastrow.h \
    WaveFunctions/partlyrestricted.h \
    WaveFunctions/slaterdeterminant.h \
    Optimization/barzilaiborwein.h \
    Optimization/gradientdescentmomentum.h

