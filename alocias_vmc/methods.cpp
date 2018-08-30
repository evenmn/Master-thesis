#include "methods.h"
#include <cmath>

void Methods::printProgressBar(std::string &flip, const float 
        &percentage, const int width, const std::string &extra) {
    /* print progressbar with percentage as output */
    flip.append(width * percentage, '-');
    flip.append(width * (1-percentage) - 1, ' ');
    std::cout << flip + "] " + std::to_string((int)(round(percentage
                    * 100.0))) + " % " + extra + "\n";
    std::cout.flush();
} // end function printProgressBar

std::string Methods::stringPos(const unsigned int &rank, const int &displ) {
    /* return a position for progress based on rank and displacement*/
    return "\033[" + std::to_string(rank + displ) + "H\r";
} // end function progressPositions

int Methods::divider(const unsigned int &index, const unsigned int &maxValue,
        const int &divisor) {
    /* divide a range in even whole number sections, return a given value
     * depending on index */
    return static_cast<int>(index < (maxValue%divisor) ?
            ceil((double)maxValue/divisor) : floor((double)maxValue/divisor));
} // end function divider

Eigen::VectorXd Methods::conjugateGradient(const Eigen::MatrixXd &A, const
        Eigen::VectorXd &rhs, const Eigen::VectorXd &x0, const unsigned int
        MAX) {
    /* solve a linear system, Ax=b with Conjugate Gradient method, return only
     * conjugate direction if number of iterations equals 1 */
    Eigen::VectorXd res = rhs - A*x0;
    if (MAX == 1) {
        return res.squaredNorm() / (res.transpose()*A*res) * res;
    } else if (MAX > 1) {
        Eigen::VectorXd p = res;
        Eigen::VectorXd xold = x0;
        Eigen::VectorXd xnew;
        Eigen::VectorXd rnew = Eigen::VectorXd::Zero(x0.size());
        double C = res.squaredNorm();
        double rInner;
        rnew << 1, 1;
        unsigned int iter = 0;
        while (iter < MAX) {
            rInner = res.squaredNorm();
            C = rInner / (p.transpose()*A*p);
            xnew = xold + C*p;
            rnew = res - C*A*p;
            if (rnew.norm() < 1e-14) {
                break;
            } // end if
            p = rnew + rnew.squaredNorm()/rInner * p;
            res = rnew;
            xold = xnew;
            iter++;
        } // end while
        return xnew;
    } // end ifeif
    return res;
} // end function conjugateGradient

double Methods::conjugateStep(const Eigen::MatrixXd &A, const
        Eigen::VectorXd &rhs, const Eigen::VectorXd &x0, const unsigned int
        MAX) {
    /* solve a linear system, Ax=b with Conjugate Gradient method and return
     * step-size */
    Eigen::VectorXd res = rhs - A*x0;
    Eigen::VectorXd p = res;
    Eigen::VectorXd rnew = Eigen::VectorXd::Zero(x0.size());
    double C = res.squaredNorm();
    double rInner;
    rnew << 1, 1;
    for (unsigned int iter = 0; iter < MAX; ++iter) {
        rInner = res.squaredNorm();
        C = rInner / (p.transpose()*A*p);
        rnew = res - C*A*p;
        if (rnew.norm() < 1e-14) {
            break;
        } // end if
        p = rnew + rnew.squaredNorm()/rInner * p;
        res = rnew;
    } // end for
    return C;
} // end function conjugateGradient

void Methods::setDisplacement(int recvc[], int displ[], const int sizes[], int
        P) {
    /* function for setting receive count and displacement used in MPI gatherv
     * and scatterv functions */
    for (int p = 0; p < P; ++p) {
        recvc[p] = sizes[p];
        displ[p] = sizes[p] * p;
    } // end forp
} // end function setDisplacement
