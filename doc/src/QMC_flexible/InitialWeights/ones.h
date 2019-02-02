#pragma once
#include "initialweights.h"

class Ones : public InitialState {
public:
    Ones(System* system, int numberOfElements);
    void setupInitialState();
};
