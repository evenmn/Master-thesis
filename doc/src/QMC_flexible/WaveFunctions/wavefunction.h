#pragma once
#include <Eigen/Dense>
#include <iostream>

class WaveFunction {

public:
    WaveFunction(class System *system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    Eigen::VectorXd getParameters() { return m_parameters; }
    virtual void updateArrays(Eigen::VectorXd positions, int pRand) = 0;
    virtual void resetArrays() = 0;
    virtual void initializeArrays(Eigen::VectorXd positions) = 0;
    virtual void updateParameters(Eigen::MatrixXd parameters) = 0;
    virtual double evaluate(Eigen::VectorXd positions) = 0;
    virtual double evaluateSqrd(Eigen::VectorXd positions) = 0;
    virtual double computeFirstDerivative(const Eigen::VectorXd positions, int k) = 0;
    virtual double computeSecondDerivative() = 0;
    virtual Eigen::VectorXd computeFirstEnergyDerivative(int k) = 0;
    virtual Eigen::VectorXd computeSecondEnergyDerivative() = 0;

    virtual ~WaveFunction() = 0;

protected:
    int     m_numberOfParticles  = 0;
    int     m_numberOfDimensions = 0;
    int     m_numberOfFreeDimensions = 0;
    int     m_numberOfParameters = 0;
    int     m_maxNumberOfParametersPerElement = 0;
    Eigen::MatrixXd m_parameters;
    //Eigen::VectorXd m_radialVector;
    //Eigen::VectorXd m_positions;
    //Eigen::VectorXd m_oldPositions;
    class System* m_system = nullptr;
};

