#include "cartesian.h"

Cartesian::Cartesian() {
} // end function constructor

Cartesian::~Cartesian() {
} // end function deconstructor

void Cartesian::setup(const unsigned int cut, const unsigned int dim) {
    s = 1;
    ms = Eigen::VectorXi(2);
    ms(0) = -1;
    ms(1) = 1;

    m_dim = dim;

    if (m_dim == 1) {
        M = Eigen::VectorXi::Zero(cut/2);
        n = Eigen::VectorXi::Zero(M.size());
        E = Eigen::VectorXi::Zero(M.size());
        states = Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>(cut, 5);
        for (unsigned int i = 0; i < cut/2; ++i) {
            M(i) = 2;
            n(i) = i;
            E(i) = i+1;
        }
        for (unsigned int i = 1; i < M.size(); ++i) {
            M(i) += M(i-1);
        } // end fori
        for (unsigned int i = 0; i < cut/2; ++i) {
            states.row(2*i)(0) = &(n(i));
            states.row(2*i)(1) = &s;
            states.row(2*i)(2) = &(ms(0));
            states.row(2*i)(3) = &(E(i));
            states.row(2*i)(4) = &(M(i));
            states.row(2*i+1)(0) = &(n(i));
            states.row(2*i+1)(1) = &s;
            states.row(2*i+1)(2) = &(ms(1));
            states.row(2*i+1)(3) = &(E(i));
            states.row(2*i+1)(4) = &(M(i));
        } // end fori
        return;
    } // end if

    // find magic numbers
    unsigned int k = 1;
    unsigned int numStates = 2;
    while (numStates < cut) {
        numStates += calculateDegeneracy(k);
        k++;
    } // end while

    M = Eigen::VectorXi::Zero(k);
    M(0) = 2;
    for (unsigned int i = 1; i < M.size(); ++i) {
        M(i) = calculateDegeneracy(i);
    } // end fori

    // set possible values for nx, ny, nz and energy E
    n = Eigen::VectorXi::Zero(M.size());
    E = Eigen::VectorXi::Zero(M.size());
    for (unsigned int i = 0; i < n.size(); ++i) {
        n(i) = i;
        E(i) = i + 1;
    } // end fori

    // allocate and set states {nx, ny, nz, s, ms, E, M}
    if (m_dim == 2) {
        states = Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>(numStates, 6);
    } else {
        states = Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>(numStates, 7);
    } // end ifelse
    unsigned int i = 0;
    unsigned int e = 0;
    Eigen::MatrixXi tmpn;
    while (i < numStates) {
        /* loop over states for energy level e */
        tmpn = Eigen::MatrixXi::Zero(M(e)/2, m_dim);
        findPrincipal(e, tmpn);
        for (unsigned int j = 0; j < tmpn.rows(); ++j) {
            if (m_dim == 2) {
                addState(i, tmpn.row(j)(0), tmpn.row(j)(1), e, 0);
                addState(i, tmpn.row(j)(0), tmpn.row(j)(1), e, 1);
            } else {
                addState(i, tmpn.row(j)(0), tmpn.row(j)(1), tmpn.row(j)(2), e, 0);
                addState(i, tmpn.row(j)(0), tmpn.row(j)(1), tmpn.row(j)(2), e, 1);
            } // end ifelse
        } // end forj
        e++;
    } // end while

    for (unsigned int i = 1; i < M.size(); ++i) {
        M(i) += M(i-1);
    } // end fori

    // set angular momenta
    angularMomenta = Eigen::VectorXi::Zero(states.rows());
    sumn();
} // end function setup

const Eigen::Ref<const Eigen::VectorXi> Cartesian::getSumn() const {
    /* return sum of n-values for all states */
    return angularMomenta;
} // end function getsumn

const int &Cartesian::getSumn(const unsigned int &i) const {
    /* return sum of n-values for state i */
    return angularMomenta(i);
} // end function getsumn

const Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>
&Cartesian::getStates() const {
    /* return states */
    return states;
} // end function getStates

const Eigen::Matrix<int*, Eigen::Dynamic, 1> Cartesian::getStates(const
        unsigned int &i) const {
    /* return state i */
    return states.row(i);
} // end function getStates

const int& Cartesian::getn(const unsigned int& i, const unsigned int& d) const
{
    /* get n in state i for dimension d */
    return *(states(i,d));
} // end function getn

const Eigen::VectorXi &Cartesian::getn() const {
    /* return vector of n-values */
    return n;
} // end function getn

const int &Cartesian::getn(const unsigned int &i) const {
    /* return value i in n-vector */
    return n(i);
} // end function getn

const Eigen::VectorXi &Cartesian::getE() const {
    /* return vector of energies */
    return E;
} // end function getE

const int &Cartesian::getE(const unsigned int &i) const {
    /* return energy of level i */
    return E(i);
} // end function getE

const Eigen::VectorXi &Cartesian::getMagic() const {
    /* return vector of magic numbers */
    return M;
} // end function getE

const int &Cartesian::getMagic(const unsigned int &i) const {
    /* return magic number for level i */
    return M(i);
} // end function getE

void Cartesian::restructureStates() {
    /* put spin-states in increasing order */
    Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> dStates =
        Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>(states.rows()/2,
                states.cols());
    Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> uStates =
        Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>(states.rows()/2,
                states.cols());
    for (unsigned int i = 0; i < states.rows(); i+=2) {
        dStates.row(i/2) = states.row(i);
        uStates.row(i/2) = states.row(i+1);
    } // end fori
    states.block(0, 0, dStates.rows(), dStates.cols()) = dStates;
    states.block(uStates.rows(), 0, uStates.rows(), uStates.cols()) = uStates;

    // recalculate angular momenta
    sumn();
} // end function restructureStates

void Cartesian::sumn() {
    /* calcualte angular momenta for each state */
    for (unsigned int i = 0; i < states.rows(); ++i) {
        angularMomenta(i) =
            Methods::refSum<int>(states.row(i).segment(0,m_dim));
    } // end fori
} // end function sumn

unsigned int Cartesian::calculateDegeneracy(const unsigned int &i) {
    /* calculate degeneracy for level i */
    return (m_dim==2 ? 2*(i+1) : (i+1)*(i+2));
} // end function calculateDegeneracy

void Cartesian::findPrincipal(const unsigned int &i, Eigen::MatrixXi &p) {
    /* find possible values for nx, ny and nz given energy level E */
    unsigned int h = 0;
    for (int j = 0; j < E(i); ++j) {
        for (int k = 0; k < E(i); ++k) {
            if (m_dim == 3) {
                for (int l = 0; l < E(i); ++l) {
                    if (j+k+l == E(i)-1) {
                        p.row(h)(0) = j;
                        p.row(h)(1) = k;
                        p.row(h)(2) = l;
                        h++;
                    } // end if
                } // end forl
            } else {
                if (j+k == E(i)-1) {
                    p.row(h)(0) = j;
                    p.row(h)(1) = k;
                    h++;
                } // end if
            } // end ifelse
        } // end fork
    } // end forj
} // end function findprincipal

void Cartesian::addState(unsigned int &i, const unsigned int j, const unsigned
        int k, const unsigned int l, const unsigned int e, const unsigned int
        ud) {
    /* create state for given level i with nx(index j, ny(index k) and nz(index
     * l) and assign to states matrix */
    states.row(i)(0) = &(n(j));
    states.row(i)(1) = &(n(k));
    states.row(i)(2) = &(n(l));
    states.row(i)(3) = &s;
    states.row(i)(4) = &(ms(ud));
    states.row(i)(5) = &(E(e));
    states.row(i)(6) = &(M(e));
    i++;
} // end function findPossiblePrincipal

void Cartesian::addState(unsigned int &i, const unsigned int j, const unsigned
        int k, const unsigned int e, const unsigned int ud) {
    /* create state for given level i with nx(index j) and ny(index k) and
     * assign to states matrix */
    states.row(i)(0) = &(n(j));
    states.row(i)(1) = &(n(k));
    states.row(i)(2) = &s;
    states.row(i)(3) = &(ms(ud));
    states.row(i)(4) = &(E(e));
    states.row(i)(5) = &(M(e));
    i++;
} // end function findPossiblePrincipal

void Cartesian::addState(unsigned int &i, const unsigned int j, const unsigned
        int e, const unsigned int ud) {
    /* create state for given level i with nx (index j) */
    states.row(i)(0) = &(n(j));
    states.row(i)(2) = &s;
    states.row(i)(3) = &(ms(ud));
    states.row(i)(4) = &(E(e));
    states.row(i)(5) = &(M(e));
    i++;
} // end function findPossiblePrincipal

void Cartesian::printStates() {
    /* print states */
    for (unsigned int i = 0; i < states.rows(); ++i) {
        std::cout << "(";
        for (int j = 0; j < states.cols(); ++j) {
            std::cout << *states.row(i)(j) << (j==states.cols()-1 ? ")" : ",");
        } // end forj
        std::cout << std::endl;
    } // end fori
    std::cout << "Number of states: " << states.rows() << std::endl;
} // end function print state
