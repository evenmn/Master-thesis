#pragma once
#include <vector>
#include <Eigen/Dense>

class WaveFunction {

public:
    WaveFunction(std::vector<class System*> system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(Eigen::MatrixXd particles) = 0;
    virtual double computeFirstDerivative(Eigen::MatrixXd particles, int k) = 0;
    virtual double computeSecondDerivative(Eigen::MatrixXd particles) = 0;
    virtual double computeFirstEnergyDerivative(Eigen::MatrixXd particles) = 0;
    virtual double computeSecondEnergyDerivative(Eigen::MatrixXd particles) = 0;

    double  TotalEvaluation(Eigen::MatrixXd particles);

    virtual ~WaveFunction() = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    std::vector<class System*> m_system;
};

