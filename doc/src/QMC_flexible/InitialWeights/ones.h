#pragma once
#include "initialweights.h"

class Ones : public InitialWeights {
public:
    Ones(System* system);
    void setupInitialWeights();
};
