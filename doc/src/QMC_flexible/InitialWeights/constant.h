#pragma once
#include "initialweights.h"

class Constant : public InitialWeights {
public:
    Constant(System* system, double factor);
    void setupInitialWeights();

private:
    double m_factor = 1;
};
