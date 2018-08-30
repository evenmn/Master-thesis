#ifndef CARTESIAN_H
#define CARTESIAN_H

#include <Eigen/Dense>

#include "../methods.h"

class Cartesian {
    private:
        int s;
        unsigned int m_dim;
        Eigen::VectorXi n, ms, E, M, m;

        Eigen::VectorXi angularMomenta;

        Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> states;

        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int);
        void findPrincipal(const unsigned int& , Eigen::MatrixXi&);
        void sumn();

        unsigned int calculateDegeneracy(const unsigned int&);
    
    public:
        Cartesian ();
        virtual ~Cartesian ();

        void setup(const unsigned int, const unsigned int);

        const Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> &getStates()
            const;
        const Eigen::Matrix<int*, Eigen::Dynamic, 1> getStates(const unsigned
                int&) const;
        const int& getn(const unsigned int&, const unsigned int&) const;
        const Eigen::VectorXi &getn() const;
        const int &getn(const unsigned int&) const;
        const Eigen::VectorXi &getE() const;
        const int &getE(const unsigned int&) const;
        const Eigen::VectorXi &getMagic() const;
        const int &getMagic(const unsigned int&) const;
        const int &getSumn(const unsigned int&) const;
        const Eigen::Ref<const Eigen::VectorXi> getSumn() const;

        void restructureStates();

        void printStates();
};

#endif /* CARTESIAN_H */
