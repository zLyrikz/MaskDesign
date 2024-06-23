// CppNumericalSolver
// //// Copyright https://github.com/PatWie/CppNumericalSolvers, MIT license
#ifndef ARMIJO_H_
#define ARMIJO_H_

#include "../meta.h"

namespace cppoptlib {

template<typename ProblemType, int Ord>
class Armijo {
public:
    using Scalar = typename ProblemType::Scalar;
    using TVector = typename ProblemType::TVector;
    /**
     * @brief use Armijo Rule for (weak) Wolfe conditiions
     * @details [long description]
     *
     * @param searchDir search direction for next update step
     * @param objFunc handle to problem
     *
     * @return step-width
     */
    static Scalar linesearch(const TVector &x, const TVector &searchDir, ProblemType &objFunc, const std::function<Scalar(const TVector&)>& _constraint) {
        const Scalar c = 1e-4;// armijo index
        const Scalar wolfe = 0.99;
        const Scalar rho = 0.9;//decrese
        Scalar alpha = 1.0;// step length
        Scalar f = 0.0;
        TVector grad(x.rows());
        const Scalar f_in = objFunc.valueAndGradient(x, grad);
        const Scalar Cache_wolfe = grad.dot(searchDir);
        const Scalar Cache = c * Cache_wolfe;
        int max_iteration = 70;//45
        int iteration = 0;
        while(iteration <= max_iteration)
        {
            if (_constraint(x + alpha * searchDir))
            {
                f = objFunc.valueAndGradient(x + alpha * searchDir, grad);
                if (f <= f_in + alpha * Cache)
                {
                    // armijo condition is met
                    // check Wolfe condition
                    const Scalar temp = grad.dot(searchDir);

                    if (temp >= wolfe * Cache_wolfe)
                    {
                        break;
                    }
                }
            }
            alpha *= rho;
            iteration++;
        }

        //if (iteration > max_iteration)
        //{
        //    cout << "BACKTRAKING iteration reach max" << endl;
        //}

        return alpha;
    }

};

template<typename ProblemType>
class Armijo<ProblemType, 2> {

 public:
    using typename ProblemType::Scalar;
    using typename ProblemType::TVector;
    using typename ProblemType::THessian;
    /**
     * @brief use Armijo Rule for (weak) Wolfe conditiions
     * @details [long description]
     *
     * @param searchDir search direction for next update step
     * @param objFunc handle to problem
     *
     * @return step-width
     */
    static Scalar linesearch(const TVector &x, const TVector &searchDir, ProblemType &objFunc) {
        const Scalar c = 0.2;
        const Scalar rho = 0.9;
        Scalar alpha = 1.0;

        Scalar f = objFunc.value(x + alpha * searchDir);
        const Scalar f_in = objFunc.value(x);
        const THessian  hessian(x.rows(), x.rows());
        objFunc.hessian(x, hessian);
        TVector grad(x.rows());
        objFunc.gradient(x, grad);
        const Scalar Cache = c * grad.dot(searchDir) + 0.5 * c*c * searchDir.transpose() * (hessian * searchDir);

        while(f > f_in + alpha * Cache) {
            alpha *= rho;
            f = objFunc.value(x + alpha * searchDir);
        }
        return alpha;
    }

};

}

#endif /* ARMIJO_H_ */
