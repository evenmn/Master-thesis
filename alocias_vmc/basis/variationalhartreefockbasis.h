#ifndef VARIATIONALHARTREEFOCKBASIS_H
#define VARIATIONALHARTREEFOCKBASIS_H

#include "../basisfunctions/cartesian.h"

#include <Eigen/Dense>

class VariationalHartreeFockBasis : public Cartesian {
    private:
        unsigned int m_dim;

    public:
        VariationalHartreeFockBasis ();
        VariationalHartreeFockBasis (unsigned int, unsigned int);
        virtual ~VariationalHartreeFockBasis ();

        void setup(unsigned int, unsigned int);
};

#endif /* VARIATIONALHARTREEFOCKBASIS_H */
