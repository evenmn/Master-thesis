#pragma once
#include "initialstate.h"

class RandomNormal : public InitialState {
public:
    RandomNormal(System* system);
    void setupInitialState();
};
