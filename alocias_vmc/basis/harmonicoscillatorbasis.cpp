#include "harmonicoscillatorbasis.h"
#include <iostream>

QuantumdotBasis::QuantumdotBasis() : Cartesian::Cartesian() {
    /* default constructor */
} // end constructor

QuantumdotBasis::QuantumdotBasis(const unsigned int cut, const unsigned int
        dim) : Cartesian::Cartesian() {
    /* initialize states */
    setup(cut, dim);
} // end constructor

QuantumdotBasis::~QuantumdotBasis() {
} // end deconstructor

void QuantumdotBasis::setup(const unsigned int cut, const unsigned int dim) {
    /* stetup basis up to a cutoff cut for dimensions dim */
    Cartesian::setup(cut, dim);
    restructureStates();
} // end function setup

void QuantumdotBasis::printDescription() {
    /* print description of printout */
    if (m_dim == 1) {
        std::cout << "(nx,s,ms,E,N)" << std::endl;
    } else if(m_dim == 2) {
        std::cout << "(nx,ny,s,ms,E,N)" << std::endl;
    } else {
        std::cout << "(nx,ny,nz,s,ms,E,N)" << std::endl;
    } // end if
} // end function

void QuantumdotBasis::printStates() {
    /* print states */
    printDescription();
    Cartesian::printStates();
    printDescription();
} // end function print state
