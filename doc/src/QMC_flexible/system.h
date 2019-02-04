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
    int                     getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    double                  getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                  getFrequency()               { return m_omega; }
    double                  getWidth()                   { return m_sigma; }
    double                  getLearningRate()            { return m_eta; }
    bool                    getInteraction()             { return m_interaction; }
    Eigen::MatrixXd         getParticles()               { return m_particles; }
    Eigen::MatrixXd         getWeights()                 { return m_parameters; }
    Eigen::MatrixXd         getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::VectorXd         getRadialVector()            { return m_radialVector; }

    double evaluateWaveFunction     (Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateWaveFunctionSqrd (Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double getKineticEnergy         ();
    void getGradient                (class WaveFunction* waveFunction, Eigen::VectorXd &TotalGradients);
    void updateParameters           (Eigen::MatrixXd &gradients);
    void calculateRadialVector      (Eigen::MatrixXd particles, Eigen::VectorXd &radialVector);
    void calculateDistanceMatrix    (Eigen::MatrixXd particles, Eigen::MatrixXd &distanceMatrix);
    std::string generate_filename   (std::string name, std::string extension);

private:
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfMetropolisSteps   = 0;
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
    Eigen::MatrixXd                     m_particles;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
};

