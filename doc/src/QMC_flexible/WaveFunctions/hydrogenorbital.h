#pragma once
#include "wavefunction.h"

class HydrogenOrbital : public WaveFunction {
public:
    HydrogenOrbital(class System* system, int elementNumber);
    void updateArrays(Eigen::VectorXd positions, int pRand);
    void resetArrays();
    void initializeArrays(Eigen::VectorXd positions);
    void updateParameters(Eigen::MatrixXd parameters);
    double evaluate();
    double evaluateSqrd();
    double computeFirstDerivative(Eigen::VectorXd positions, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    double calculateRadialVectorElement(Eigen::VectorXd particles, int par);
    Eigen::VectorXd calculateRadialVector(Eigen::VectorXd particles);

private:
    double m_alpha = 0;
    int m_elementNumber = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_oldRadialVector;
};
