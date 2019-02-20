#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, int elementNumber);
    void updateArrays(Eigen::VectorXd positions, int pRand);
    void resetArrays();
    void initializeArrays(Eigen::VectorXd positions);
    void updateParameters(Eigen::MatrixXd parameters);
    double evaluate(Eigen::VectorXd positions);
    double evaluateSqrd(Eigen::VectorXd positions);
    double computeFirstDerivative(const Eigen::VectorXd positions, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    Eigen::MatrixXd calculateDistanceMatrix(Eigen::VectorXd particles);
    double calculateDistanceMatrixElement(int i, int j, Eigen::VectorXd particles);
    Eigen::MatrixXd calculateDistanceMatrixCross(int par, Eigen::VectorXd particles);

    //double f(int i, int j);
    Eigen::MatrixXd f();
    double g(int i, int j, int k, int l);
    double beta(int i, int j);
    double gamma();

private:
    int m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_oldDistanceMatrix;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
};
