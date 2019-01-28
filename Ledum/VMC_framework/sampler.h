#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    int     m_acceptenceRatio = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_dE = 0;
    double  m_dEE = 0;
    double  m_SqrdE = 0;
    class System* m_system = nullptr;
};
