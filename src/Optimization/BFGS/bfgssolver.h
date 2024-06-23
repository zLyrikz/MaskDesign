// CppNumericalSolver
// //// Copyright https://github.com/PatWie/CppNumericalSolvers, MIT license
#include <iostream>
#include <fstream>
#include <Eigen/LU>
#include "isolver.h"
#include "../linesearch/armijo.h"

#ifndef BFGSSOLVER_H_
#define BFGSSOLVER_H_

namespace cppoptlib {

template<typename ProblemType>
class BfgsSolver : public ISolver<ProblemType, 1> {
  public:
    using Superclass = ISolver<ProblemType, 1>;
    using typename Superclass::Scalar;
    using typename Superclass::TVector;
    using typename Superclass::THessian;

    void minimize(ProblemType &objFunc, TVector & x0) {
        const size_t DIM = x0.rows();
        THessian H = THessian::Identity(DIM, DIM);
        TVector grad(DIM);
        TVector x_old = x0;
        TVector x_initial = x0;
        this->m_current.reset();
        objFunc.gradient(x0, grad);
        TVector x_best = x0;
        Scalar f_best = objFunc.value(x0);


        do {
            TVector searchDir = -1 * H * grad;
            // check "positive definite"
            Scalar phi = grad.dot(searchDir);

            // positive definit ?
            if ((phi > 0) || (phi != phi)) {
                //cout << "****************** hessian not possitive, reset the hessian approximation *******************************" << endl;
                // no, we reset the hessian approximation
                H = THessian::Identity(DIM, DIM);
                searchDir = -1 * grad;
            }

            const Scalar rate = Armijo<ProblemType, 1>::linesearch(x0, searchDir, objFunc, constraint_) ;

            x0 = x0 + rate * searchDir;
            if (!constraint_(x0))
            {
                // not feasible
                // restart
                //cout << endl << "###################### not feasible, restart #################################" << endl << endl;
                x0 = x_best;
                objFunc.gradient(x0, grad);
                H = THessian::Identity(DIM, DIM);
                continue;
            }
            else 
            {
                Scalar current_f = objFunc.value(x0);
                if (current_f < f_best)
                {
                    x_best = x0;
                    f_best = current_f;
                }
            }

            TVector grad_old = grad;
            objFunc.gradient(x0, grad);
            TVector s = rate * searchDir;
            TVector y = grad - grad_old;

            const Scalar rho = 1.0 / y.dot(s);
            H = H - rho * (s * (y.transpose() * H) + (H * y) * s.transpose()) + rho * (rho * y.dot(H * y) + 1.0)
                * (s * s.transpose());


            //if( (x_old-x0).template lpNorm<Eigen::Infinity>() < 1e-7  )
            //    break;
            x_old = x0;



            ++this->m_current.iterations;
            this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
            this->m_status = checkConvergence(this->m_stop, this->m_current);

            //std::cout << "iter: " << this->m_current.iterations/*<< " f = " <<  objFunc.value(x0)*/ << " ||g||_inf " << this->m_current.gradNorm <<" step length="<< s.norm() << std::endl;
            //std::cout << "   x-x_initial = " << std::endl;

        } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));

        x0 = x_best;

    }

    std::function<Scalar(const TVector&)> constraint_;

};

}
/* namespace cppoptlib */

#endif /* BFGSSOLVER_H_ */
