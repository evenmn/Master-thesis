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
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (std::vector<class WaveFunction*> waveFunctionVector);
    void setInitialState            (class InitialState* initialState);
    void setOptimizer               (class Optimization* optimization);
    double evaluateWaveFunction     (const Eigen::MatrixXd &particles);
    double getKineticEnergy         (const Eigen::MatrixXd &particles);
    double getGradient              (const Eigen::MatrixXd &particles);
    class WaveFunction*                          getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*                           getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                               getSampler()        { return m_sampler; }
    class Optimization*                          getOptimizer()      { return m_optimizer; }
    Eigen::MatrixXd                              getParticles()      { return m_particles; }
    int    getNumberOfParticles()       { return m_numberOfParticles; }
    int    getNumberOfDimensions()      { return m_numberOfDimensions; }
    int    getNumberOfMetropolisSteps() { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    class WaveFunction*             m_waveFunction = nullptr;
    std::vector<class WaveFunction*> m_waveFunctionVector;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    class Optimization*             m_optimizer = nullptr;
    Eigen::MatrixXd                 m_particles; // = Eigen::MatrixXd();
};

