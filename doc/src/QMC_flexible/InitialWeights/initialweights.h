#pragma once
#include <Eigen/Dense>

class InitialWeights {
public:
    InitialWeights(class System* system);
    virtual void setupInitialWeights() = 0;
    Eigen::MatrixXd getWeights() { return m_parameters; }

    virtual ~InitialWeights() = 0;

protected:
    class System* m_system = nullptr;
    int              m_numberOfDimensions    = 0;
    int              m_numberOfParticles     = 0;
    int              m_numberOfElements      = 0;
    int              m_maxNumberOfParametersPerElement = 0;
    Eigen::MatrixXd  m_parameters;
};
