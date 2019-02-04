#pragma once
#include <Eigen/Dense>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep, int stepNumber);
    void printOutputToTerminal();
    void computeAverages(Eigen::MatrixXd &gradients);
    double getEnergy()          { return m_energy; }
    void getEnergyGradient(double EL_avg, Eigen::MatrixXd grad_tot, Eigen::MatrixXd gradE_tot, Eigen::MatrixXd &gradients);

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfParticles = 0;
    int     m_numberOfDimensions = 0;
    int     m_stepNumber = 0;
    int     m_acceptenceRatio = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    Eigen::MatrixXd  m_dE;
    Eigen::MatrixXd  m_dEE;
    double  m_SqrdE = 0;
    class System* m_system = nullptr;
};
