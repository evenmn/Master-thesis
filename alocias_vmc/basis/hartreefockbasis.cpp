#include "hartreefockbasis.h"

HartreeFockBasis::HartreeFockBasis() : Cartesian() {
    /* default constructor */
} // end constructor

HartreeFockBasis::HartreeFockBasis(unsigned int cut, unsigned int dim) :
    Cartesian() {
    /* default constructor */
    setup(cut, dim);
} // end constructor

HartreeFockBasis::~HartreeFockBasis() {
} // end deconstructor

void HartreeFockBasis::setup(unsigned int cut , unsigned int dim) {
    /* initiate states */
    Cartesian::setup(cut, dim);
    Cartesian::restructureStates();
} // end function setup
