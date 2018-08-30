#ifndef HARTREEFOCKBASIS_H
#define HARTREEFOCKBASIS_H

#include "../basisfunctions/cartesian.h"

#include <Eigen/Dense>

class HartreeFockBasis : public Cartesian {
    private:
        unsigned int m_dim;

    public:
        HartreeFockBasis ();
        HartreeFockBasis (unsigned int, unsigned int);
        virtual ~HartreeFockBasis ();
        
        void setup(unsigned int, unsigned int);
};

#endif /* HARTREEFOCKBASIS_H */
