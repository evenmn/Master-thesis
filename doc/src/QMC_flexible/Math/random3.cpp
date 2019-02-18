#include "random3.h"
#include <iostream>

namespace mt {
  #include "mersenne-twister/mersenne-twister.h"
}

double Random3::nextGaussian(double mean, double standardDeviation) {
    return 1;
}

int Random3::nextInt(int upperLimit) {
    return upperLimit * int(mt::rand_u32());
}

double Random3::nextDouble() {
    return mt::rand_u32();
}
