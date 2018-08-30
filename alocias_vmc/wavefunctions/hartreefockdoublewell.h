#ifndef HARTREEFOCKDOUBLEWELL_H
#define HARTREEFOCKDOUBLEWELL_H

#include "../basis/hartreefockdoublewellbasis.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Slater;

class HartreeFockDoubleWell : public HartreeFockDoubleWellBasis {
    private:
        bool m_interaction, isFull;

        unsigned int m_numBasis;

        double omega, omegaSq, sqrtOmega;
        
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
        void setHermiteNormalizations();
        void setHartreeFockBasisNormalizations();

    public:
        HartreeFockDoubleWell(Slater*);
        virtual ~HartreeFockDoubleWell();
        
        void setParameters(const Eigen::VectorXd&);
        void initializeParameters(const double& , Eigen::MatrixXd);
        void setInteraction(bool);
        
        std::string setupDone();

        double potentialEnergy();

        double calculateWavefunction(const unsigned int&, const unsigned int&);
        
        double gradientExpression(const unsigned int&, const int&, const
                unsigned int&);
        const Eigen::VectorXd& laplacianExpression(const unsigned int&, const
                unsigned int&);
        double variationalDerivativeExpression(const unsigned int&, const
                unsigned int&, const double&);

        void reSetAll();
        void initializeMatrices();
        void set(const Eigen::MatrixXd& newPositions);
        void update(const Eigen::VectorXd&,  const unsigned int&);
        void reset(const unsigned int&);
        void acceptState(const unsigned int&);

    protected:
        std::vector<double(HartreeFockDoubleWell::*)(const unsigned int&, const unsigned int&,
                const double&)> variationalDerivativeFunctionList;
        
};

#endif /* HARTREEFOCKDOUBLEWELL_H */
