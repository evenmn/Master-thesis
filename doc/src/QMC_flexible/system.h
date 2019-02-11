#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

class System {
public:
    bool metropolisStep             ();
    void runMetropolisSteps         (int numberOfMetropolisSteps, int numberOfIterations);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setNumberOfFreeDimensions  ();
    void setMaxNumberOfParametersPerElement (int maxNumberOfParametersPerElement);
    void setNumberOfWaveFunctionElements (int numberOfWaveFunctionElements);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setFrequency               (double omega);
    void setWidth                   (double sigma);
    void setLearningRate            (double eta);
    void setInteraction             (bool interaction);
    void setGradients               ();
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setInitialState            (class InitialState* initialState);
    void setInitialWeights          (class InitialWeights* initialWeights);
    //void setOptimizer               (class Optimization* optimizer);
    void setWaveFunction            (std::vector<class WaveFunction*> waveFunctionVector);

    class WaveFunction*     getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*      getHamiltonian()             { return m_hamiltonian; }
    class Sampler*          getSampler()                 { return m_sampler; }
    //class Optimization*     getOptimizer()               { return m_optimizer; }
    int                     getNumberOfParticles()       { return m_numberOfParticles; }
    int                     getNumberOfDimensions()      { return m_numberOfDimensions; }
    int                     getNumberOfFreeDimensions()  { return m_numberOfFreeDimensions; }
    int                     getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    int                     getMaxNumberOfParametersPerElement() { return m_maxNumberOfParametersPerElement; }
    int                     getNumberOfWaveFunctionElements() { return m_numberOfWaveFunctionElements; }
    double                  getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                  getFrequency()               { return m_omega; }
    double                  getWidth()                   { return m_sigma; }
    double                  getLearningRate()            { return m_eta; }
    bool                    getInteraction()             { return m_interaction; }
    Eigen::VectorXd         getParticles()               { return m_particles; }
    Eigen::VectorXd         getPositions()               { return m_positions; }
    Eigen::MatrixXd         getWeights()                 { return m_parameters; }
    Eigen::MatrixXd         getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::VectorXd         getRadialVector()            { return m_radialVector; }

    double evaluateWaveFunction     (Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateWaveFunctionSqrd (Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double getKineticEnergy         ();
    Eigen::VectorXd getGradient(class WaveFunction* waveFunction);
    Eigen::MatrixXd updateParameters();
    Eigen::VectorXd calculateRadialVector      (Eigen::VectorXd particles);
    Eigen::MatrixXd calculateDistanceMatrix    (Eigen::VectorXd particles);
    std::string generate_filename   (std::string name, std::string extension);

private:
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfFreeDimensions    = 0;
    int                                 m_numberOfMetropolisSteps   = 0;
    int                                 m_numberOfWaveFunctionElements = 0;
    int                                 m_maxNumberOfParametersPerElement = 0;
    bool                                m_interaction               = false;
    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    double                              m_omega                     = 1.0;
    double                              m_sigma                     = 1.0;
    double                              m_eta                       = 0.1;
    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    //class Optimization*                 m_optimizer                 = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
    Eigen::VectorXd                     m_particles;
    Eigen::VectorXd                     m_positions;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
};

