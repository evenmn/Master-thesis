#ifndef MTLS_H 
#define MTLS_H 

#include <eigen3/Eigen/Dense>
#include <memory>
#include <iostream>

class MTLS {
    private:
        MTLS () {};
        virtual ~MTLS () {};

        template<typename T> static inline void updateBrackets(T& al, const T&
                ac, T& au, T& fl, T& fc, T& fu, T& dl, T& dc, T& du) {
            /* update end-points in bracket in More-Thuente linesearch */
            if (fc > fl) {
                /* squeeze end-point from upper-side, changing the end-points
                 * of the bracket to ac and au */
                au = ac, fu = fc, du = dc;
            } else {
                /* squeeze end-point from lower-side */
                if (dc*(al-ac) > 0) {
                    /* change the end-points of the bracket to al and ac */
                    al = ac, fl = fc, dl = dc;
                } else {
                    /* flip (interchange) end-points and set ac to best */
                    au = al, fu = fl, du = dl;
                    al = ac, fl = fc, dl = dc;
                } // end ifelse
            } // end ifelse
        } // end function updateBrackets

        template<typename T> static inline void updateBrackets(T& al, const T&
                ac, T& au, const T& mfl, const T& mfc, const T& mdc, T& fl, T&
                fc, T& fu, T& dl, T& dc, T& du) {
            /* update end-points in bracket in More-Thuente linesearch with
             * auxiliary bounds */
            if (mfc > mfl) {
                /* squeeze end-point from upper-side, changing the end-points
                 * of the bracket to ac and au */
                au = ac, fu = fc, du = dc;
            } else {
                /* squeeze end-point from lower-side */
                if (mdc*(al-ac) > 0) {
                    /* change the end-points of the bracket to al and ac */
                    al = ac, fl = fc, dl = dc;
                } else {
                    /* flip (interchange) end-points and set ac to best */
                    au = al, fu = fl, du = dl;
                    al = ac, fl = fc, dl = dc;
                } // end ifelse
            } // end ifelse
        } // end function updateBrackets

        template<typename T> static inline T stepMoreThuente(const T& al, const
                T& ac, const T& au, const T& fl, const T& fc, const T& fu,
                const T& dl, const T& dc, const T& du, const T& aMin, const T&
                aMax, bool& isBracketed, bool& bound) {
            /* find new trial step in More-Thuente linesearch */
            T aNew = ac;

            // find new trial step
            if (fc > fl) {
                /* case 1: function value higher */
                T cubicInt = cubicPolMin(al, ac, fl, fc, dl, dc);
                T quadflfcdl = quadPolf0f1g0Min(al, ac, fl, fc, dl);
                if (fabs(cubicInt - al) < fabs(quadflfcdl - al)) {
                    /* choose point closest to al (best estimate so-far) */
                    aNew = cubicInt;
                } else {
                    /* choose midpoint as a compromise */
                    aNew =  0.5 * (cubicInt + quadflfcdl);
                } // end ifelse

                // ac is within brackets and curvature suggests cubic is within
                // bracket bounds
                isBracketed = true;
                bound = true;
            } else if (dc*dl < 0) {
                /* case 2: function value lower and point between al and ac */
                T cubicInt = cubicPolMin(al, ac, fl, fc, dl, dc);
                T quadfldldc = quadPolf0g0g1Min(al, ac, dl, dc);
                if (fabs(cubicInt - ac) >= fabs(quadfldldc - ac)) {
                    /* choose point closest to al (smallest estimate so-far) */
                    aNew = cubicInt;
                } else {
                    /* case: quadratic is better than cubic */
                    aNew = quadfldldc;
                } // end ifelse

                // ac is within brackets, but curvature suggests cubic is
                // outside of bounds 
                isBracketed = true;
                bound = false;
            } else if (fabs(dc) < fabs(dl)) {
                /* case 3: function value lower, but cubic interpolator may not
                 * have a minimum */
                T cubicInt = cubicPolMin(al, ac, fl, fc, dl, dc, aMin, aMax);
                T quadfldldc = quadPolf0g0g1Min(al, ac, dl, dc);
                if (fabs(cubicInt - ac) < fabs(quadfldldc - ac)) {
                    /* choose cubic minimum if it is within brackets */
                    aNew = cubicInt; 
                } else {
                    aNew = quadfldldc;
                } // end ifelse

                // no information about brackets, but cubic is bound
                bound = true;
            } else {
                /* case 4: function decreasing rapidly in search direction, but
                 * interpolators do not give any more information */
                if (isBracketed) {
                    aNew = cubicPolMin(ac, au, fc, fu, dc, du);
                } else if (al < ac) {
                    aNew = aMax;
                } else {
                    aNew = aMin;
                } // end ifeifelse

                // no info on brackets, cubic is not bound
                bound = false;
            } // end ifeifeifeif

            // make sure trial step is lower than maximum step
            if (aNew > aMax) {
                aNew = aMax;
            } // end if
            
            // make sure trial step is higher than minimum step
            if (aNew < aMin) {
                aNew = aMin;
            } // end if

            return aNew;
        } // end function stepMoreThuente
        
        template<typename T> static inline int sgn(T x) {
            /* general sign function */
            return (x > T(0)) - (x < T(0));
        } // end function sgn
        
        template<typename T> static inline T quadPolf0f1g0Min(T x_0, T x_1, T
                f_0, T f_1, T g_0) {
            T diff_x = x_0 - x_1;
            T diff_f = f_0 - f_1;
            return -0.5*(2*diff_f*x_0 - diff_x*g_0*(x_0 + x_1))/(-diff_f +
                    diff_x*g_0);
        } // end function quadPolf0f1g0Min

        template<typename T> static inline T quadPolf0f1g1Min(T x_0, T x_1, T
                f_0, T f_1, T g_1) {
            T diff_x = x_0 - x_1;
            T diff_f = f_0 - f_1;
            return -0.5*(-2*diff_f*x_1 + diff_x*g_1*(x_0 + x_1))/(diff_f -
                    diff_x*g_1);
        } // end function quadPolf0f1g1Min

        template<typename T> static inline T quadPolf0g0g1Min(T x_0, T x_1, T
                g_0, T g_1) {
            return -1.0*(-g_0*x_1 + g_1*x_0)/(g_0 - g_1);
        } // end function quadPolf0g0g1Min

        template<typename T> static inline T quadPolf1g0g1Min(T x_0, T x_1, T
                g_0, T g_1) {
            return -1.0*(-g_0*x_1 + g_1*x_0)/(g_0 - g_1);
        } // end function quadPolf1g0g1Min

        template<typename T> static inline T cubicPolMin(T x0, T x1, T f0, T
                f1, T g0, T g1) {
            /* expression for minimum of cubic interpolator given in Nocedal &
             * Wright Numerical Optimization, equation 3.59, page 59. */
            T fdiff = f0 - f1;
            T xdiff = x1 - x0;
            T d1 = g0 + g1 + 3*(fdiff / xdiff);
            T d2 = sgn(xdiff) * sqrt(d1*d1 - g1*g0);
            return x1 - xdiff * ((g1 + d2 - d1) / (g1 - g0 + 2*d2));
        } // end function cubicPolMin

        template<typename T> static inline T cubicPolMin(T x0, T x1, T f0, T
                f1, T g0, T g1, T aMin, T aMax) {
            /* expression for minimum of cubic interpolator given in Nocedal &
             * Wright Numerical Optimization, equation 3.59, page 59. With the
             * added test for if cubic tends to infinity in the direction of
             * current minimmizer x1 (that is g1) */
            T fdiff = f0 - f1;
            T xdiff = x1 - x0;
            T d1 = g0 + g1 + 3*(fdiff / xdiff);
            T d2 = sgn(xdiff) * sqrt(d1*d1 - g1*g0);
            T ratio = (g1 + d2 - d1) / (g1 - g0 + 2*d2);
//             T cubicInt =  x1 - xdiff * ((g1 + d2 - d1) / (g1 - g0 + 2*d2));
            if ((ratio < 0) && (d2 != 0)) {
                /* choose cubic if it tends to infinity in the direction of
                 * minimizer */
                return x1 - xdiff * ratio;
            } else if (d1 < 0) {
                /* return maximum if cubic points opposite to minimizer and is
                 * beyond maximum step */
                return aMin;
            } else {
                /* return minimum if cubic is points opposite to minimizer and
                 * is beyond minimum step */
                return aMax;
            } // end ifeifelse
        } // end function cubicPolMin

        template<typename T> static inline T min(const T& a1, const T& a2) {
            /* return minimum og a1 and a2 */
            return (a1 < a2 ? a1 : a2);
        } // end function min
        
        template<typename T> static inline T max(const T& a1, const T& a2) {
            /* return maximum og a1 and a2 */
            return (a1 > a2 ? a1 : a2);
        } // end function max

        template<typename T, typename F, typename G, typename... Args> class Dummy {
            /* dummy class for containing calculation function func and
             * derivative derFunc */
            public:
                Dummy() {};
                virtual ~Dummy() {};

                Eigen::Matrix<T, Eigen::Dynamic, 1> derivative;

                F f;
                G df;

                T eva(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, Args...
                        ags) {
                    derivative = df(x);
                    return f(x, ags...);
                } // end function eva

                const Eigen::Matrix<T, Eigen::Dynamic, 1>& der() const {
                    return derivative;
                } // end function der
        }; // end class Dummy
        
        template<typename T, typename F, typename G, typename... Args> class
            DummyRaw {
            /* dummy class for containing calculation function func and
             * derivative derFunc that return plain arrays */
            private:
                unsigned int m_xSize;

            public:
                DummyRaw(const unsigned int& xSize) {
                    m_xSize = xSize;
                };
                virtual ~DummyRaw() {};

                Eigen::Matrix<T, Eigen::Dynamic, 1> derivative;

                F f;
                G df;

                T eva(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, Args...
                        ags) { 
                    derivative = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic,
                               1>>(df(x.data()), m_xSize, 1);
                    return f(x.data(), ags...);
                } // end function eva

                const Eigen::Matrix<T, Eigen::Dynamic, 1>& der() const {
                    return derivative;
                } // end function der
        }; // end class DummyRaw

    public:
        template<typename T, typename F, typename G, typename... Args> static
            inline T linesearchMoreThuente(T* searchDirection, T* x0, const T
                    f0, F func, G derFunc, Args... args) {
                /* override in case raw arrays are given and an object for
                 * calculation function and parameters struct are both not
                 * given */
            unsigned int xSize = sizeof(x0) / sizeof(T);
            std::shared_ptr<DummyRaw<T,F,G,Args...>> d =
                std::make_unique<DummyRaw<T,F,G,Args...>>(xSize);

            Eigen::Matrix<T, Eigen::Dynamic, 1> searchDirE =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic,
                1>>(searchDirection, xSize, 1);
            Eigen::Matrix<T, Eigen::Dynamic, 1> x0E =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(x0, xSize, 1);

            d->f = func;
            d->df = derFunc;
            d->derivative = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic,
                1>>(derFunc(x0), xSize, 1);

            return linesearchMoreThuente(searchDirE, x0E, f0, d.get(),
                    &DummyRaw<T,F,G,Args...>::eva,
                    &DummyRaw<T,F,G,Args...>::der, args...);
        } // end function lineseearchMoreThuente
        
        template<typename P, typename T, typename F, typename G, typename...
            Args> static inline T linesearchMoreThuente(P* inputs, T*
                    searchDirection, T* x0, const T f0, F func, G derFunc,
                    Args... args) {
            /* override in case raw arrays are given and an object for
             * calculation function is not given */

            // find size of array and create dummy object
            size_t xSize = sizeof(x0) / sizeof(T);
            std::shared_ptr<DummyRaw<T,F,G,Args...>> d =
                std::make_unique<DummyRaw<T,F,G,Args...>>(xSize);

            // convert plain arrays to eigen vectors
            Eigen::Matrix<T, Eigen::Dynamic, 1> searchDirE =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic,
                1>>(searchDirection, xSize, 1);
            Eigen::Matrix<T, Eigen::Dynamic, 1> x0E =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(x0, xSize, 1);

            // set calculation and derivative function
            d->f = func;
            d->df = derFunc;
            d->derivative = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic,
                1>>(derFunc(x0), xSize, 1);

            return linesearchMoreThuente(inputs, searchDirE, x0E, f0, d.get(),
                    &DummyRaw<T,F,G,Args...>::eva,
                    &DummyRaw<T,F,G,Args...>::der, args...);
        } // end function lineseearchMoreThuente

        template<typename T, typename F, typename G, typename...  Args>
            static inline T linesearchMoreThuente(const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& searchDirection, const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& x0, const T f0, F func, G derFunc,
                    Args... args) {
            /* override in case an object for calculation function and
             * parameters struct are both not given */
            
            // create object and set functions inside dummy (essentially using
            // Dummy as wrapper)
            std::shared_ptr<Dummy<T,F,G,Args...>> d =
                std::make_unique<Dummy<T,F,G,Args...>>();
            d->f = func;
            d->df = derFunc;
            d->derivative = derFunc(x0);

            return linesearchMoreThuente(searchDirection, x0, f0, d.get(),
                    &Dummy<T,F,G,Args...>::eva , &Dummy<T,F,G,Args...>::der,
                    args...);
        } // end function linesearchMoreThuente
        
        template<typename P, typename T, typename F, typename G, typename
            EigenVecT=Eigen::Matrix<T, Eigen::Dynamic, 1>, typename...  Args>
            static inline T linesearchMoreThuente(P* inputs, const
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& searchDirection, const
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& x0, const T f0, F
                    func, G derFunc, Args... args) {
            /* override in case an object for calculation function is not given
             * */

            // create object and set functions inside dummy (essentially using
            // Dummy as wrapper)
            std::shared_ptr<Dummy<T,F,G,Args...>> d =
                std::make_unique<Dummy<T,F,G,Args...>>();
            d->f = func;
            d->df = derFunc;
            d->derivative = derFunc(x0);

            return linesearchMoreThuente(inputs, searchDirection, x0, f0,
                    d.get(), &Dummy<T, F, G, Args...>::eva , &Dummy<T, F, G,
                    Args...>::der, args...);
        } // end function linesearchMoreThuente
        
        template<typename T, class U, typename F, typename G, typename...
            Args> static inline T linesearchMoreThuente(const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& searchDirection, const Eigen::Matrix<T,
                    Eigen::Dynamic, 1>& x0, const T f0, const U& funcObj, F
                    func, G derFunc, Args... args) {
            /* override to set default parameters and then run linesearch in
             * case a struct of parameters is not given */
            struct Params {
                /* struct of default parameters */
                T maxIterations = 100; // maximum number of iterations
                T mu = 0.5; // step scaling factor, (0<mu<=1/2<eta)
                T eta = 1.0; // termination parameter (0<eta<1)
                T delta = 4.0; // delta: scaling of step updating ([1.1,4.0])
                T bisectWidth = 0.66; // extrapolation tolerance for bounds
                T bracketTol = 1e-14; // termination tolerance for brackets
                T aMin0 = 0.0; // lower bound for step (aMin>=0 and aMin<aMax)
                T aMax0 = 100.0; // upper bound for step (aMax>aMin)
            } defaultParams; // end struct Params

            return linesearchMoreThuente(&defaultParams, searchDirection, x0,
                    f0, funcObj, func, derFunc, args...);
        } // end function linesearchMoreThuente
        
        template<typename P, typename T, class U, typename F, typename G,
            typename ...Args> static inline T linesearchMoreThuente(P* inputs,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p, const
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& x0, const T f0, const
                    U& funcObj, F func, G derFunc, Args...  args) {
            /* perform line search as described in paper 'Line Search
             * Algorithms with Guaranteed Sufficient Decrease',
             * doi:10.1145/192115.192132, by Jorge J. More and David J.
             * Thuente.
             *
             * explanation of input arguments:
             *  inputs: struct build with variables:
             *      maxIterations: maximum number of iterations
             *      mu: step scaling factor, 0<mu<=1/2<eta (notice mu != 0)
             *      eta: termination parameter 0<eta<1 (notice noninclusivety)
             *      delta: scaling of step updating (in interval [1.1,4.0])
             *      bisectWidth: extrapolation tolerance for bounds
             *      bracketTol: termination tolerance for brackets 
             *      aMin0: lower bound for step with aMin>=0 and aMin<aMax 
             *      aMax0: lower bound for step with aMin>=0 and aMin<aMax 
             *  searchDir: search direction
             *  x0: current set of parameters
             *  f0: current function value
             *  
             *  funcObj: object containing evaluate function func
             *  func: function for calculating new function value
             *  derFunc: function for grabbing derivative
             *  args(optional): arguments to be passed to func (not derFunc) */

            // derivative with respect to step
            T d0 = p.dot((funcObj->*derFunc)());

            // products reused in calculation of auxiliary function and in
            // termination test
            T d0mu = inputs->mu*d0;
            T d0eta = fabs(d0)*inputs->eta;

            // set initial bounds
            T al = 0.0; // end-point, smallest step so-far (bracket)
            T au = 0.0; // other end-point for step (bracket)
            T ac = 1.0; // current(or trial) step

            // function- and derivative values
            T fl = f0; // best so-far
            T dl = d0; // derivative at al
            T fu = f0; // value at other end-point
            T du = d0; // derivative at other end-point
            T fc = f0; // trial function value
            T dc = d0; // trial step
          
            // current boundaries
            T aMin = inputs->aMin0;
            T aMax = inputs->aMax0;

            // width of bracket (used for forcing sufficient decrease)
            T width = aMax - aMin;
            T prevWidth = 2 * width;

            // value for determining if step is within brackets and if
            // auxiliary function is needed and if cubic interpolation point is
            // within bounds
            bool isBracketed = false;
            bool auxiliaryNeeded = true;
            bool bound = false;

            for (unsigned int k = 0; k < inputs->maxIterations; ++k) {
                /* perform iterations */

                //////////////////
                // ORDER BOUNDS //
                //////////////////
                if (isBracketed) {
                    /* order aMin and aMax if ac is within good brackets */
                    aMin = min(al, au);
                    aMax = max(al, au);
                } else {
                    /* set minimum to current smallest and safeguard maximum */
                    aMin = al;
                    aMax = ac + inputs->delta * (ac - al);
                } // end ifelse

                // make sure step is within boundaries [aMin0, aMax0]
                if (ac > inputs->aMax0) {
                    /* make sure step is lower than maximum step */
                    ac = inputs->aMax0;
                } // end if
                if (ac < inputs->aMin0) {
                    /* make sure step is higher than minimum step */
                    ac = inputs->aMin0;
                } // end if

                // take care of wierd behavior (due to either round off errors
                // or instability with interpolation points)
                if ((isBracketed && ((ac <= aMin) || (ac >= aMax))) ||
                        (isBracketed && ((aMax-aMin) <=
                                         inputs->bracketTol*aMax))) {
                    /* set step to current smallest */
                    ac = al;
                } // end if
             
                //////////////////////////
                // CALCULATE NEW VALUES //
                //////////////////////////
                fc = (funcObj->*func)(x0 + ac*p, args...);
                dc = p.dot((funcObj->*derFunc)());
//                 fc = ((*funcObj).*func)(x0 + ac*p, args...);
//                 dc = p.dot(((*funcObj).*derFunc)());

                //////////////////////////
                // TEST FOR CONVERGENCE //
                //////////////////////////
                T ftest = f0 + d0mu*ac;

                if ((ac == inputs->aMax0) && ((fc <= ftest) && (dc <= d0eta)))
                {
                    /* step is at maximum value */
                    break;
                } // end if
                
                if ((ac == inputs->aMin0) && ((fc > ftest) || (dc >= d0eta))) {
                    /* step is at minimum value */
                    break;
                } // end if
                
                if (isBracketed && ((aMax - aMin) <= (inputs->bracketTol *
                                aMax))) {
                    /* bracket is small enough */
                    break;
                } // end if
                
                if ((fc <= ftest) && (fabs(dc) <= d0eta)) {
                    /* sufficient decrease and strong curvature conditions */
                    break;
                } // end if

                if (auxiliaryNeeded && ((fc <= ftest) && (dc > 0))) {
                    /* check if auxiliary is needed any further */
                    auxiliaryNeeded = false;
                } // end if

                /////////////////////////////////////////////
                // FIND NEW TRIAL STEP AND UPDATE BRACKETS //
                /////////////////////////////////////////////
                T trialStep = ac;
                if (auxiliaryNeeded) {
                    /* use auxiliary function */
                    T mfl = fl - (f0 + al*d0mu);
                    T mdl = dl - d0mu;
                    T mfu = fu - (f0 + au*d0mu);
                    T mdu = du - d0mu;
                    T mfc = fc - (f0 + ac*d0mu);
                    T mdc = dc - d0mu;

                    trialStep = stepMoreThuente(al, ac, au, mfl, mfc, mfu, mdl,
                            mdc, mdu, aMin, aMax, isBracketed, bound);
                    updateBrackets(al, ac, au, mfl, mfc, mdc, fl, fc, fu, dl,
                            dc, du);
                } else {
                    /* use modified algorithm */
                    trialStep = stepMoreThuente(al, ac, au, fl, fc, fu, dl, dc,
                            du, aMin, aMax, isBracketed, bound);
                    updateBrackets(al, ac, au, fl, fc, fu, dl, dc, du);
                } // end ifelse
               
                // refine the trial step (stay careful with interpolation)
                if (isBracketed && bound) {
                    /* make sure point is within boundary (ac, au) and not too
                     * close to au */
                    T bVal = al + inputs->bisectWidth * (au - al);
                    if (al < au) {
                        if (bVal < trialStep) {
                            trialStep = bVal;
                        } // end if
                    } else {
                        if (trialStep < bVal) {
                            trialStep = bVal;
                        } // end if
                    } // end ifelse
                } // end if

                // update current step
                ac = trialStep;

                ///////////////////////////////
                // FORCE SUFFICIENT DECREASE //
                ///////////////////////////////
                if (isBracketed) {
                    /* force reduction of bracket size */
                    if ((inputs->bisectWidth * prevWidth) <= fabs(au - al)) {
                        /* set trial to mid-point if outside */
                        ac = 0.5 * (au + al);
                    } // end if
                    prevWidth = width;
                    width = fabs(au - al);
                } // end if
            } // end fork

            return ac;
        } // end function linesearchMoreThuente
}; // end class MTLS

#endif /* MTLS_H */
