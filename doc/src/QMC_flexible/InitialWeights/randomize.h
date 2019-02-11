#pragma once
#include "initialweights.h"

class Randomize : public InitialWeights {
public:
    Randomize(System* system);
    void setupInitialWeights();
};
