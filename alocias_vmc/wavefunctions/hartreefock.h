#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <memory>
#include <string>

#include "../basis/hartreefockbasis.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Slater;

class HartreeFock : public HartreeFockBasis {
    private:
        double omega, omegaSq, sqrtOmega;
        unsigned int m_numBasis;

        bool m_interaction, isFull;
        
        Eigen::VectorXd m_laplacianSumVec;
        
        Eigen::MatrixXd m_SnewPositions, m_SoldPositions,
            m_SWavefunctionMatrix, m_SoldWavefunctionMatrix;

        Slater* slater;
        
        Eigen::SparseMatrix<double> m_C;
        Eigen::ArrayXd m_hermiteNormalizations;
        
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            m_hermite3DMatrix, m_oldHermite3DMatrix;

        void setHermite3DMatrix(const unsigned int&);

        void setBasisWavefunction(const unsigned int&);
        void setBasisWavefunction();
        void setHermiteNormalizations(Eigen::MatrixXd&);

    public:
        HartreeFock(Slater*);
        virtual ~HartreeFock();
        
        void setParameters(const Eigen::VectorXd&);
        void initializeParameters(const double&, const unsigned int&,
                Eigen::MatrixXd);
        void setInteraction(bool);
        
        std::string setupDone();
        void checkIfFullShell();

        double potentialEnergy();

        double calculateWavefunction(const unsigned int&, const unsigned int&);
        
        double gradientExpression(const unsigned int&, const int&, const
                unsigned int&);
        const Eigen::VectorXd& laplacianExpression(const unsigned int&, const
                unsigned int&);
        double variationalDerivativeExpression(const unsigned int&, const
                unsigned int&);

        void reSetAll();
        void initializeMatrices();
        void set(const Eigen::MatrixXd& newPositions);
        void update(const Eigen::VectorXd&,  const unsigned int&);
        void reset(const unsigned int&);
        void acceptState(const unsigned int&);

    protected:
        std::vector<double(HartreeFock::*)(const unsigned int&, const unsigned
                int&, const double&)> variationalDerivativeFunctionList;
        
};

#endif /* HARTREEFOCK_H */
