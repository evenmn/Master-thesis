#ifndef METHODS_H
#define METHODS_H

#include <eigen3/Eigen/Dense>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>

class Methods {
    private:
        Methods () {};
        virtual ~Methods () {};

    public:
        static void printProgressBar(std::string&, const float&, const int,
                const std::string& extra="");
        static std::string stringPos(const unsigned int&, const int&);
        static int divider(const unsigned int&, const unsigned int&, const
                int&);
        static Eigen::VectorXd conjugateGradient(const Eigen::MatrixXd&, const
                Eigen::VectorXd&, const Eigen::VectorXd&, const unsigned
                int=100); 
        static double conjugateStep(const Eigen::MatrixXd&, const
                Eigen::VectorXd&, const Eigen::VectorXd&, const unsigned
                int=100);
        static void setDisplacement(int[], int[], const int[], int);

        template<typename T=double, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>>
        static unsigned int argMin(const EigenVectorXT& vec) {
            /* return index of coefficient of minimum value */
            unsigned int minIdx = 0;
            T minValue = vec(0);
            for (unsigned int i = 1; i < vec.size(); ++i) {
                if (vec(i) < minValue) {
                    minValue = vec(i);
                    minIdx = i;
                } // end if
            } // end fori
            return minIdx;
        } // end function argMin
        
        template<typename T=double, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>>
        static unsigned int argMax(const EigenVectorXT& vec) {
            /* return index of coefficient of maximum value */
            unsigned int maxIdx = 0;
            T maxValue = vec(0);
            for (unsigned int i = 1; i < vec.size(); ++i) {
                if (vec(i) > maxValue) {
                    maxValue = vec(i);
                    maxIdx = i;
                } // end if
            } // end fori
            return maxIdx;
        } // end function argMax

        template<typename T=double, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>>
        static bool strongWolfeCondition(
                const T& beta,
                const EigenVectorXT& newDerivative,
                const EigenVectorXT& oldDerivative,
                const EigenVectorXT& searchDirection) {
            /* check that Strong-Wolfe condition, return true if met and false
             * if not */
            if (innerProd(newDerivative, searchDirection) < (beta *
                        innerProd(oldDerivative, searchDirection))) {
                /* conditions are met */
                return true;
            } else {
                /* conditions are not met */
                return false;
            } // end ifelse
        } // end function strongWolfeCondition

        template<
            typename T=double, 
            typename EigenVectorXT=Eigen::Matrix<T, Eigen::Dynamic, 1>, 
            typename EigenMatrixXT=Eigen::Matrix<T, Eigen::Dynamic,
            Eigen::Dynamic>> 
        static T innerProd(const EigenVectorXT& vec, const EigenMatrixXT& mat)
        {
        /* calculate and return inner-product between Eigen column-vector vec
         * and matix mat */
            return vec.transpose() * mat * vec;
        } // end function innerProd
        
        template<
            typename T=double, 
            typename EigenVectorXT=Eigen::Matrix<T, Eigen::Dynamic, 1>>
        static T innerProd(const EigenVectorXT& v, const EigenVectorXT& u) {
            return v.transpose() * u;
        } // end function innerProd

        template<typename T, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>> static void steepestDescent(const T& step,
                    EigenVectorXT& parameters, const EigenVectorXT& gradient) {
            /* minimize variational parameters with steepest descent */
            parameters -= step * gradient;
        } // end function steepestDescent

        template<typename T> static inline void switcher(T& a, T& b) {
            /* switch a and b */
            T tmpa = a;
            a = b;
            b = tmpa;
        } // end function switcher

        template<typename T> static inline int sgn(T x) {
            /* general sign function */
            return (x > T(0)) - (x < T(0));
        } // end function sgn

        template<typename T> static inline void updateMatrixInverse(const T&
                Mold, const T& Mnew, const T& MoldInv, Eigen::Ref<T> MnewInv,
                const double& ratio, const unsigned int i) {
            /* update inverse of matrix when only row i has changed */
            MnewInv.setZero();
            for (unsigned int k = 0; k < MnewInv.rows(); ++k) {
                for (unsigned int j = 0; j < MnewInv.rows(); ++j) {
                    for (unsigned int l = 0; l < MnewInv.rows(); ++l) {
                        if (j==i) {
                            MnewInv(k,j) += Mold(i,l)*MoldInv(l,j);
                        } else {
                            MnewInv(k,j) -= Mnew(i,l)*MoldInv(l,j);
                        } // end ifelse
                    } // end forl
                    MnewInv(k,j) *= MoldInv(k,i)/ratio;
                    MnewInv(k,j) = ((j!=i) ? MnewInv(k,j)+MoldInv(k,j) :
                            MnewInv(k,j));
                } // end forj
            } // end fork
        } // end function updateMatrixInverse 

        template<typename T> static inline double determinantRatio(const T&
                newElement, const T& oldInverse, const unsigned int i) {
            /* Calculate ratio of determinants when only one row has changed */
            double R = 0;
            for (unsigned int j = 0; j < oldInverse.rows(); ++j) {
                R += newElement(i,j) * oldInverse(j,i);
            } // end fori
            return R;
        } // end function determinantRatio

        template<typename T, typename U, typename V> bool withinRange(const T&
                value, const U& area, const V& mid) {
            return (((value > (mid-area)) && (value < (mid+area))) ? true :
                    false);
        } // end function within range

        template<typename U> U vecDiff(const Eigen::Matrix<const U,
                Eigen::Dynamic, Eigen::Dynamic>& vector) {
            return vector(0) - vector(1);
        } // end function vecDiff
        
        template<typename U> U rowwiseVecdiff(const Eigen::Matrix<const U,
                Eigen::Dynamic, Eigen::Dynamic>& matrix) {
            Eigen::VectorXd diffVector = Eigen::VectorXd::Zero(matrix.rows());
            for (unsigned int i = 0; i < matrix.rows(); ++i) {
                /* loop over rows */
                diffVector(i) = vecdiff(matrix.row(i));
            } // end fori
        } // end function rowwiseVecdiff

        template<typename U, typename T> static inline U inverseFraction(const
                T &val) {
            /* return 1/val */
            return 1./val;
        } // end function inverseFraction

        template <typename T> static inline T refSum(const Eigen::Matrix<T*,
                Eigen::Dynamic, 1>& vec) {
            /* sum values of a vector of pointers(de-reference and sum) */
            T sum = 0;
            for (unsigned int i = 0; i < vec.size(); ++i) {
                sum += *(vec(i));
            } // end fori
            return sum;
        } // end function refSum

        template <typename T> static void breaker(T val) {
            /* print val */
            std::cout << std::setprecision(14) << val << " ";
            return;
        } // end template function print1

        struct expand {
            /* struct for expanding function pack */
            template<typename... T> expand(T&&...) {}
        };
        template<typename... Args> static void breaker(int max, Args... vals) {
            /* call breaker for all vals (poor-mans debugging to a new
             * level...) */
            static int b = 0;
            if (b >= max) {
                exit(1);
            } else {
                expand{0, (breaker(vals), 0)...};
                std::cout << std::endl;
                b++;
            } // end if
        } // end function breaker

        template<typename T> static inline void print1(T arg) {
            std::cout << arg << " ";
        } // end function print1
        
        template<typename... Args> static inline void sepPrint(Args... args) {
            /* print args with spaces */
            expand{0, (print1(args), 0)...};
            std::cout << "\n";
        } // end function sepPrint

        template<typename T> static inline unsigned int reMapIndex(const T&
                val, const unsigned int &max) {
            T inv = 1./((T)max);
            for (unsigned int i = 0; i < max; ++i) {
                /* loop over number og dimensions */
                if ((val>=i*inv) && (val<(inv*(i+1)))) {
                    /* return index (0,1,...) if value is within range */
                    return i;
                } // end ifelse
            } // end ford
            return max-1;
        } // end function reMapIndex
        
        template<typename T> static void matrix3DSetZero(T &matrix, const
                unsigned int &i, const unsigned int &d) {
            /* set all values in row i along dimension d to zero */
            for (unsigned int j = 0; j < matrix.cols(); ++j) {
                matrix(i,j)(d) = 0;
            } // end forj
        } // end function matrix3DSetZero
        
        template<typename T> static void matrix3DSetZero(T &matrix, const
                unsigned int &i) {
            /* set all values along row i in 3D matrix to zero */
            for (unsigned int j = 0; j < matrix.cols(); ++j) {
                matrix(i,j).setZero();
            } // end forj
        } // end function matrix3DSetZero
        
        template<typename T> static void matrix3DSetZero(T &matrix) {
            /* set all values in matrix to zero */
            for (int i = 0; i < matrix.rows(); ++i) {
                matrix3DSetZero(matrix, i);
            } // end fori
        } // end function matrix3DSetZero

        template <typename T> static void symSwap(T& A, T& B, const unsigned
                int& i) {
            /* swap row and column i of matrix A and B */
            A.row(i) = B.row(i);
            A.col(i).head(i) = B.col(i).head(i);
            A.col(i).tail(A.rows()-(i+1)) = B.col(i).tail(A.rows()-(i+1));
        } // end function symSwap
        
        template <typename T> static void symSwapNoDiag(T& A, T& B, const
                unsigned int& i) {
            /* swap row and column i of matrix A and B, ignoring diagonal
             * elements */
            A.row(i).head(i) = B.row(i).head(i);
            A.row(i).tail(A.rows()-(i+1)) = B.row(i).tail(A.rows()-(i+1));
            A.col(i).head(i) = B.col(i).head(i);
            A.col(i).tail(A.rows()-(i+1)) = B.col(i).tail(A.rows()-(i+1));
        } // end function symSwap
        
        template <typename T> static void symSwapNoDiag(T& A, T& B, const
                unsigned int& i, const unsigned int& span) {
             /* Override which takes precalculated span as argument */
            A.row(i).head(i) = B.row(i).head(i);
            A.row(i).tail(span) = B.row(i).tail(span);
            A.col(i).head(i) = B.col(i).head(i);
            A.col(i).tail(span) = B.col(i).tail(span);
        } // end function symSwap

        template<typename T> static inline T min(const T& a1, const T& a2) {
            /* return minimum og a1 and a2 */
            return (a1 < a2 ? a1 : a2);
        } // end function min
        
        template<typename T> static inline T max(const T& a1, const T& a2) {
            /* return maximum og a1 and a2 */
            return (a1 > a2 ? a1 : a2);
        } // end function max

        template<typename T> static inline T var(const Eigen::Ref<const
                Eigen::Array<T, Eigen::Dynamic, 1>> array) {
            /* calculate and return variance of an Eigen array */
            return (array - array.mean()).pow(2).mean();
        } // end function variance
        
        template<typename T> static inline T var(const Eigen::Ref<const
                Eigen::Array<T, Eigen::Dynamic, 1>> array, const T& mean) {
            /* calculate and return variance of an Eigen array, override with
             * precalculated mean */
            return (array - mean).pow(2).mean();
        } // end function variance
        
        template<typename T> static inline Eigen::Matrix<T, Eigen::Dynamic, 1>
            cwisePow(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, const
                    Eigen::VectorXi& pows) {
                /* raise each element in x to the power of corresponding
                 * elements in pows */
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = x;
            for (unsigned int i = 0; i < x.size(); ++i) {
                res(i) = pow(x(i), pows(i));
            } // end fori
            return res;
        } // end function cwisePow

        template<typename T, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>> static inline void BFGSInverse(Eigen::Matrix<T,
                    Eigen::Dynamic, Eigen::Dynamic> &F, const EigenVectorXT&
                    xNew, const EigenVectorXT& xOld, const EigenVectorXT&
                    fDerivativeNew, const EigenVectorXT& fDerivativeOld) {
                /* update inverse matrix F with
                 * Broyden-Fletcher-Goldfarb-Shanno method (basically applied
                 * the Sherman-Morrison formula to the original BFGS formula)
                 * */
                EigenVectorXT y = fDerivativeNew - fDerivativeOld;
                EigenVectorXT s = xNew - xOld;

                T sy = s.transpose() * y;
                if (sy < 1e-14) {
                    /* make sure parameters are not to close to each-other */
                    return;
                } // end if

                T denom = 1. / sy;
                F += ((sy + y.transpose()*F*y) * s*s.transpose())
                    * (denom*denom) - (F*y*s.transpose() + s*y.transpose()*F) *
                    denom;
        } // end function BFGSInverse
        
        template<typename T, typename EigenVectorXT=Eigen::Matrix<T,
            Eigen::Dynamic, 1>> static inline void BFGS(Eigen::Matrix<T,
                    Eigen::Dynamic, Eigen::Dynamic>& F, const EigenVectorXT&
                    xNew, const EigenVectorXT& xOld, const EigenVectorXT&
                    fDerivativeNew, const EigenVectorXT& fDerivativeOld) {
            /* update F with the Broyden-Fletcher-Goldfarb-Shanno method */
            EigenVectorXT xDiff = xNew - xOld;
            if (xDiff.norm() < 1e-14) {
                /* make sure parameters are stable */
                return;
            } // end if

            EigenVectorXT fDerivativeDiff = fDerivativeNew - fDerivativeOld;
            if (fDerivativeDiff.norm() < 1e-14) {
                /* make sure derivatives are stable */
                return;
            } // end if

            EigenVectorXT FdotXdiff = F * xDiff;

            F += fDerivativeDiff*fDerivativeDiff.transpose() /
                fDerivativeDiff.dot(xDiff) - FdotXdiff*FdotXdiff.transpose() /
                (xDiff.dot(FdotXdiff));
        } // end function BFGS

        template<typename T, class U, typename F, typename... Args> static
            inline T backtrackingLinesearch(T step, const T& stepChange, const
                    T& tol, const T& oldValue, const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& xold, const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& searchDirection, const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& derivative, const U& funcObj, F func,
                    Args... args) {
                /* perform backtracking line search with function func in
                 * funcObject along search direction 'searchDirection', return
                 * new set of parameters updated with step as satisfied by
                 * Armijo-Goldstein condition. */

                // perform linesearch
                Eigen::Matrix<T, Eigen::Dynamic, 1> xnew;
                T tolerance = tol * searchDirection.dot(derivative);
                do {
                    /* perform backtracking line-search until step satisfies
                     * Armijo-Goldstein condition */
                    step *= stepChange;
                    xnew = xold + searchDirection*step;
                } while (!((funcObj->*func)(xnew, args...) <= (oldValue +
                                step*tolerance))); // end do while

                return step;
        } // end function backtrackingLinesearch
}; // end class Methods

#endif /* METHODS_H */
