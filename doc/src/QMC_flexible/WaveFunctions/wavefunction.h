#pragma once
#include <Eigen/Dense>
#include <iostream>

class WaveFunction {

public:
    WaveFunction(class System *system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    Eigen::VectorXd getParameters() { return m_parameters; }
    virtual double evaluate(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) = 0;
    virtual double evaluateSqrd(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) = 0;
    virtual double computeFirstDerivative(int k) = 0;
    virtual double computeSecondDerivative() = 0;
    virtual void computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) = 0;
    virtual void computeSecondEnergyDerivative(Eigen::VectorXd &gradients) = 0;

    virtual ~WaveFunction() = 0;

protected:
    int     m_numberOfParticles  = 0;
    int     m_numberOfDimensions = 0;
    int     m_numberOfParameters = 0;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_particles;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    class System* m_system = nullptr;
};

