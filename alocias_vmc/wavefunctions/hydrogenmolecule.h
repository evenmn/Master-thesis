#ifndef HYDROGENMOLECULE_H
#define HYDROGENMOLECULE_H

#include <Eigen/Dense>

class Hydrogenmolecule {
    private:
        unsigned int dim, m_numParticles;
        double alpha, beta, R;

        Eigen::MatrixXd positionsProtons, m_newPositions, m_oldPositions,
            m_newDistances, m_oldDistances;
        
        Eigen::MatrixXd m_gradientMatrix;
        Eigen::MatrixXd m_jastrowGradientMatrix;

        Eigen::VectorXd m_parameters, m_firstDerivativesParameters;

        Eigen::MatrixXd m_newWavefunctionMatrix, m_oldWavefunctionMatrix;

        void setDistances();
        void setDistances(const unsigned int&);

        double calculateWavefunction(const unsigned int);
        double calculateWavefunction(const unsigned int, const unsigned int);
    
    public:
        Hydrogenmolecule (double, const unsigned int);
        virtual ~Hydrogenmolecule ();

        double wavefunctionRatio();
        void gradient();
        void gradient(const unsigned int&);
        void gradient(const unsigned int&, const unsigned int&);
        double laplacian();

        double jastrow(const unsigned int&);
        void calculateJastrowGradient();
        void calculateJastrowGradient(const unsigned int&);
        double jastrowLaplacian();

        void update(const double&, const unsigned int&, const unsigned int&);
        void update(const Eigen::VectorXd&, const unsigned int&);
        void reset(const unsigned int&, const unsigned int&);
        void reset(const unsigned int&);
      
        void resetGradient(const unsigned int&, const unsigned int&);
        void resetGradient(const unsigned int&);
        void resetGradient();
        
        void resetJastrowGradient(const unsigned int&, const unsigned int&);
        void resetJastrowGradient(const unsigned int&);
        void resetJastrowGradient();

        void acceptState();
        void acceptState(const unsigned int&);
        void acceptState(const unsigned int&, const unsigned int&);
        
        void acceptGradient();
        void acceptGradient(const unsigned int&);
        void acceptGradient(const unsigned int&, const unsigned int&);
        
        void acceptJastrowGradient();
        void acceptJastrowGradient(const unsigned int&);
        void acceptJastrowGradient(const unsigned int&, const unsigned int&);

        void set(const Eigen::MatrixXd&);
        void rset(const Eigen::MatrixXd&, unsigned int);
        
        void setParameters(const Eigen::VectorXd&);
        void setVariationalDerivatives(bool);
        void setNumParticles(const unsigned int);

        void initialize(const Eigen::MatrixXd&);

        void setWavefunction();
        void setWavefunction(const unsigned int&);

        double potentialEnergy();
        double kineticEnergy(bool);

        const unsigned int &getDimension() const;
        const double &getNewPosition(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getNewPosition(const
                unsigned int&) const;
        const Eigen::MatrixXd &getNewPositions() const;
        const double &getOldPosition(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getOldPosition(const
                unsigned int&) const;
        const Eigen::MatrixXd &getOldPositions() const ;
        const Eigen::VectorXd &getDerivativesParameters() const;
        const double &getGradient(const unsigned int&, const unsigned int&)
            const;
        const Eigen::VectorXd getGradient(const unsigned int&) const;
        const Eigen::MatrixXd &getGradientMatrix() const;
        const double &getJastrowGradient(const unsigned int&, const unsigned
                int&) const;
        const Eigen::VectorXd getJastrowGradient(const unsigned int&) const;
        const Eigen::MatrixXd &getJastrowGradientMatrix() const;
        const unsigned int &getNumberOfParticles() const;
        const Eigen::VectorXd &getVariationalDerivatives() const;
};

typedef Hydrogenmolecule Wavefunction;

#endif /* HYDROGENMOLECULE_H */
