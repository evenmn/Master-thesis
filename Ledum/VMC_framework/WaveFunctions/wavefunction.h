#pragma once
#include <vector>
#include <Eigen/Dense>

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(Eigen::MatrixXd particles) = 0;
    virtual double computeDerivative(Eigen::MatrixXd particles) = 0;
    virtual double computeEnergyDerivative(Eigen::MatrixXd particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

