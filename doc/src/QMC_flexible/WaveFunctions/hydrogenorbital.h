#pragma once
#include "wavefunction.h"

class HydrogenOrbital : public WaveFunction {
public:
    HydrogenOrbital(class System* system, double beta);
    void updateArrays(Eigen::VectorXd positions, int pRand);
    void resetArrays();
    void initializeArrays(Eigen::VectorXd positions);
    void updateParameters(Eigen::MatrixXd parameters);
    double evaluate();
    double computeDerivative(Eigen::MatrixXd positions);
    double computeEnergyDerivative(Eigen::MatrixXd positions);

    Eigen::VectorXd calculateRadialVector(Eigen::VectorXd particles);

private:
    double m_beta = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_oldRadialVector;
};
