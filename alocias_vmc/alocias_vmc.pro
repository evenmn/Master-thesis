TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    slaterjastrow.cpp \
    slater.cpp \
    resampler.cpp \
    methods.cpp \
    wavefunctions/variationalhartreefockdoublewell.cpp \
    wavefunctions/variationalhartreefock.cpp \
    wavefunctions/testwavefunction.cpp \
    wavefunctions/hydrogenmolecule.cpp \
    wavefunctions/hartreefockdoublewell.cpp \
    wavefunctions/hartreefock.cpp \
    wavefunctions/harmonicoscillator.cpp \
    jastrow/rbmjastrow.cpp \
    jastrow/padenqs.cpp \
    jastrow/padejastrow.cpp \
    jastrow/expnqs.cpp \
    examples/run.cpp \
    basisfunctions/cartesian.cpp \
    basis/variationalhartreefockdoublewellbasis.cpp \
    basis/variationalhartreefockbasis.cpp \
    basis/testwavefunctionbasis.cpp \
    basis/testwavebasis.cpp \
    basis/hartreefockdoublewellbasis.cpp \
    basis/hartreefockbasis.cpp \
    basis/harmonicoscillatorbasis.cpp

HEADERS += \
    vmc.h \
    slaterjastrow.h \
    slater.h \
    resampler.h \
    MTLS.h \
    minimizer.h \
    methods.h \
    importanceSampling.h \
    hasmemfunc.h \
    bruteforce.h \
    accumulate.h \
    wavefunctions/variationalhartreefockdoublewell.h \
    wavefunctions/variationalhartreefock.h \
    wavefunctions/testwavefunction.h \
    wavefunctions/hydrogenmolecule.h \
    wavefunctions/hartreefockdoublewell.h \
    wavefunctions/hartreefock.h \
    wavefunctions/harmonicoscillator.h \
    jastrow/rbmjastrow.h \
    jastrow/padenqs.h \
    jastrow/padejastrow.h \
    jastrow/expnqs.h \
    hermite/hermite.h \
    hermite/gausscontracted.h \
    examples/run.h \
    basisfunctions/dwc.h \
    basisfunctions/cartesian.h \
    basis/variationalhartreefockdoublewellbasis.h \
    basis/variationalhartreefockbasis.h \
    basis/testwavefunctionbasis.h \
    basis/testwavebasis.h \
    basis/hartreefockdoublewellbasis.h \
    basis/hartreefockbasis.h \
    basis/harmonicoscillatorbasis.h \
    yaml-cpp/include/yaml-cpp/yaml.h \
    yaml-cpp/include/yaml-cpp/traits.h \
    yaml-cpp/include/yaml-cpp/stlemitter.h \
    yaml-cpp/include/yaml-cpp/parser.h \
    yaml-cpp/include/yaml-cpp/ostream_wrapper.h \
    yaml-cpp/include/yaml-cpp/null.h \
    yaml-cpp/include/yaml-cpp/noncopyable.h \
    yaml-cpp/include/yaml-cpp/mark.h \
    yaml-cpp/include/yaml-cpp/exceptions.h \
    yaml-cpp/include/yaml-cpp/eventhandler.h \
    yaml-cpp/include/yaml-cpp/emitterstyle.h \
    yaml-cpp/include/yaml-cpp/emittermanip.h \
    yaml-cpp/include/yaml-cpp/emitterdef.h \
    yaml-cpp/include/yaml-cpp/emitter.h \
    yaml-cpp/include/yaml-cpp/emitfromevents.h \
    yaml-cpp/include/yaml-cpp/dll.h \
    yaml-cpp/include/yaml-cpp/binary.h \
    yaml-cpp/include/yaml-cpp/anchor.h \
    yaml-cpp/binary.h \
    yaml-cpp/anchor.h \
    yaml-cpp/yaml.h \
    yaml-cpp/traits.h \
    yaml-cpp/stlemitter.h \
    yaml-cpp/parser.h \
    yaml-cpp/ostream_wrapper.h \
    yaml-cpp/null.h \
    yaml-cpp/noncopyable.h \
    yaml-cpp/mark.h \
    yaml-cpp/exceptions.h \
    yaml-cpp/eventhandler.h \
    yaml-cpp/emitterstyle.h \
    yaml-cpp/emittermanip.h \
    yaml-cpp/emitterdef.h \
    yaml-cpp/emitter.h \
    yaml-cpp/emitfromevents.h \
    yaml-cpp/dll.h

# MPI Settings
release {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
