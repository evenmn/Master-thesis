#pragma once

class Random2 {
public:
    static int    nextInt(int upperLimit);
    static double nextDouble();
    static double nextGaussian(double mean, double standardDeviation);
};
