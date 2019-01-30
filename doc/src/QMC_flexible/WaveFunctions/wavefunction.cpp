#include "wavefunction.h"
#include <iostream>

WaveFunction::WaveFunction(std::vector<System *> system) {
    m_system = system;
}

WaveFunction::~WaveFunction() {};
