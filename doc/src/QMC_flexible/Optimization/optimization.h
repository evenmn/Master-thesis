#pragma once
#include <Eigen/Dense>
#include <vector>

class Optimization {
public:
    Optimization(class System* system);
    virtual Eigen::MatrixXd updateParameters() = 0;
    virtual Eigen::MatrixXd getAllImmediateGradients() = 0;
    virtual ~Optimization() = 0;

protected:
    class System* m_system = nullptr;
    Eigen::VectorXd m_positions;
    int m_numberOfFreeDimensions = 0;
    int m_numberOfStepsAfterEquilibrium = 0;
    int m_numberOfWaveFunctionElements = 0;
    int m_maxNumberOfParametersPerElement = 0;
    double m_eta = 1;
    std::vector<class WaveFunction*> m_waveFunctionVector;
};
