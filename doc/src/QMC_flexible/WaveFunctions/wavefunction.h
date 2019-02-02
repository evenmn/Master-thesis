#pragma once
#include <Eigen/Dense>
#include <iostream>

class WaveFunction {

public:
    WaveFunction(class System *system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    Eigen::VectorXd getParameters() { return m_parameters; }
    virtual double evaluate(Eigen::MatrixXd particles) = 0;
    virtual double computeFirstDerivative(int k) = 0;
    virtual double computeSecondDerivative() = 0;
    virtual double computeFirstEnergyDerivative() = 0;
    virtual double computeSecondEnergyDerivative() = 0;

    virtual ~WaveFunction() = 0;

protected:
    int     m_numberOfParticles  = 0;
    int     m_numberOfDimensions = 0;
    int     m_numberOfParameters = 0;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_particles;
    class System* m_system = nullptr;
};

