#pragma once
#include "initialweights.h"

class Ones : public InitialWeights {
public:
    Ones(System* system, int numberOfElements);
    void setupInitialWeights();
};
