#ifndef VARIATIONALHARTREEFOCKDOUBLEWELL_H
#define VARIATIONALHARTREEFOCKDOUBLEWELL_H

#include "../basis/variationalhartreefockdoublewellbasis.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Slater;

class VariationalHartreeFockDoubleWell : 
    public VariationalHartreeFockDoubleWellBasis {
    private:
        bool m_interaction;
        
        unsigned int m_numBasis;

        double omega, omegaSq, sqrtOmega, aw, sqrtaw;
        
        Eigen::VectorXd m_laplacianSumVec;

        Slater* slater;
        
        Eigen::MatrixXd m_SnewPositions, m_SoldPositions,
            m_SWavefunctionMatrix, m_SoldWavefunctionMatrix;
        
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
        VariationalHartreeFockDoubleWell(Slater*);
        virtual ~VariationalHartreeFockDoubleWell();
        
        void setParameters(const Eigen::VectorXd&);
        void initializeParameters(const double& w, Eigen::MatrixXd);
        void setInteraction(bool);

        double potentialEnergy();
        
        std::string setupDone();

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
        std::vector<double(VariationalHartreeFockDoubleWell::*)(const unsigned
                int&, const unsigned int&, const double&)>
            variationalDerivativeFunctionList;
        
};

#endif /* VARIATIONALHARTREEFOCKDOUBLEWELL_H */
