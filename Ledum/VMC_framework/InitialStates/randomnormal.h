#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState {
public:
    RandomNormal(System* system, int numberOfDimensions, int numberOfParticles);
    void setupInitialState();
};
