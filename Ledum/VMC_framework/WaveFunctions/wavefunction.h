#pragma once
#include <vector>
#include <Eigen/Dense>

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    Eigen::VectorXd getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    //int     m_numberOfDimensions = 0;
    //double  m_alpha = 0;
    Eigen::VectorXd m_parameters = Eigen::VectorXd();
    class System* m_system = nullptr;
};

