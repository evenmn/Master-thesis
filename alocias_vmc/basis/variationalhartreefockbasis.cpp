#include "variationalhartreefockbasis.h"

VariationalHartreeFockBasis::VariationalHartreeFockBasis() : Cartesian() {
    /* default constructor */
} // end constructor

VariationalHartreeFockBasis::VariationalHartreeFockBasis(unsigned int cut,
        unsigned int dim) : Cartesian() {
    /* default constructor */
    setup(cut, dim);
} // end constructor

VariationalHartreeFockBasis::~VariationalHartreeFockBasis() {
} // end deconstructor

void VariationalHartreeFockBasis::setup(unsigned int cut , unsigned int dim) {
    /* initiate states */
    Cartesian::setup(cut, dim);
    Cartesian::restructureStates();
} // end function setup
