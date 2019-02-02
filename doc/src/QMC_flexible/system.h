#pragma once
#include <Eigen/Dense>
#include <vector>

class System {
public:
    bool metropolisStep             ();
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setInteraction             (bool interaction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setInitialState            (class InitialState* initialState);
    void setInitialWeights          (class InitialWeights* initialWeights);
    void setOptimizer               (class Optimization* optimization);
    void setWaveFunction            (std::vector<class WaveFunction*> waveFunctionVector);

    class WaveFunction*     getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*      getHamiltonian()             { return m_hamiltonian; }
    class Sampler*          getSampler()                 { return m_sampler; }
    class Optimization*     getOptimizer()               { return m_optimizer; }
    Eigen::MatrixXd         getParticles()               { return m_particles; }
    Eigen::MatrixXd         getWeights()                 { return m_parameters; }
    int                     getNumberOfParticles()       { return m_numberOfParticles; }
    int                     getNumberOfDimensions()      { return m_numberOfDimensions; }
    int                     getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    double                  getEquilibrationFraction()   { return m_equilibrationFraction; }
    bool                    getInteraction()             { return m_interaction; }

    double evaluateWaveFunction     (Eigen::MatrixXd particles);
    double getKineticEnergy         ();
    double getGradient              ();

private:
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfMetropolisSteps   = 0;
    bool                                m_interaction               = false;
    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    class Optimization*                 m_optimizer                 = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
    Eigen::MatrixXd                     m_particles;
    Eigen::MatrixXd                     m_parameters;
};

