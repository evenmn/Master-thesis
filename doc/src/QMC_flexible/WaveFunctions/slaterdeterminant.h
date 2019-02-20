#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system, int elementNumber);
    Eigen::MatrixXd list();
    double updateElement(Eigen::VectorXd positions, double basis(double, int), int i, int j);
    Eigen::VectorXd updateRow(Eigen::VectorXd positions, double basis(double, int), int i);
    Eigen::MatrixXd updateMatrix(Eigen::VectorXd positions, double basis(double, int));
    Eigen::VectorXd dA_row(Eigen::VectorXd positions, int k);
    Eigen::MatrixXd dA_matrix(Eigen::VectorXd positions);

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
private:
    int     m_elementNumber     = 2;
    double  m_omega             = 0;
    double  m_alpha             = 0;
    int m_numberOfOrbitals      = 0;
    int m_numberOfParticlesHalf = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::MatrixXd m_D_up;
    Eigen::MatrixXd m_D_dn;
    Eigen::MatrixXd m_D_upOld;
    Eigen::MatrixXd m_D_dnOld;
    Eigen::MatrixXd m_dD_up;
    Eigen::MatrixXd m_dD_dn;
    Eigen::MatrixXd m_dD_upOld;
    Eigen::MatrixXd m_dD_dnOld;
};
