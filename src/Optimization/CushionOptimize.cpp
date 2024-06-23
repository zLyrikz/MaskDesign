#include "CushionOptimize.h"
#include "ConstraintFunctor.h"

#include "ObjectiveForCppoptlib.h"
#include "../Optimization/BFGS/bfgssolver.h"
#include "../FEM/WrapperAlglibFunction.h"
#include "../Utility/TypeConvert.h"
#include "../Alglib/optimization.h"
#include "meta.h"
#include <iostream>
#include <regex>
using std::cout;
using std::endl;
CushionOptimize::CushionOptimize(CrossSectionOptimizeInfo* _cross_section_info)
{
	cross_section_info_ = _cross_section_info;
}

CushionOptimize::~CushionOptimize()
{
}

void CushionOptimize::BfgsOptimize(CushionSurface& _optimal_cushion, int _max_iteration, double _epsilon_gradient)
{
    // initialize the objective
    ObjectiveForCppoptlib objective;
    objective.setCrossSectionInfo(cross_section_info_);

    // Set initial guess.
    int dim_x = cross_section_info_->GetDimX();
    const double* initial_x = cross_section_info_->GetInitialX();
    Eigen::VectorXd x_best(dim_x);
    for (int iX = 0; iX < dim_x; ++iX)
    {
        x_best[iX] = initial_x[iX];
    }

    // choose a solver
    cppoptlib::BfgsSolver<ObjectiveForCppoptlib> solver;
    cppoptlib::Criteria<double> stop_critria;
    stop_critria.gradNorm = _epsilon_gradient;
    stop_critria.iterations = _max_iteration;
    solver.setStopCriteria(stop_critria);
    ConstraintFunctor constraint(cross_section_info_);
    solver.constraint_ = constraint;
    //std::cout << "f initial " << objective(x_best) << std::endl;

    // and minimize the function
    solver.minimize(objective, x_best);

    //// print argmin
    //std::cout << "argmin      " << x_best.transpose() << std::endl;
    //std::cout << "f in argmin " << objective(x_best) << std::endl;
    cross_section_info_->ChangeCushionFromX(x_best.data(), _optimal_cushion, cross_section_info_->GetSampleV());

}

bool CushionOptimize::GetFEMVisualizationInfo(
    std::vector<double>& _area_distribution_color, std::vector<double>& _force_distribution_color, std::vector<double>& _sigmoid_color,
    TetrahedralMesh& _tetrahedral_cushion, Mesh& _deformed_cushion_surface,
    Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale) const
{
    const LinearTetrahedralFemEvaluation* fem_evaluator = cross_section_info_->GetFemEvaluator();
    if (fem_evaluator->IsSimulationGood())
    {
        fem_evaluator->AreaDistributionColor(_area_distribution_color);
        fem_evaluator->ForceDistributionColor(_force_distribution_color);
        fem_evaluator->SigmoidColor(_sigmoid_color);

        fem_evaluator->GetDeformedCushion(_tetrahedral_cushion, _deformed_cushion_surface);

        _tetrahedral_cushion.FindBoundaryVertex();
        fem_evaluator->GetTetrahedralCushionForces(_tetrahedral_cushion, _force_mesh, _node_force_magnitude, _scale);

        //std::vector<double> cushion_surface_nodes_force_magnitude;
        //fem_evaluator->GetFreeSurfaceNodesForce(cushion_surface_nodes_force_magnitude);
        return true;
    }
    else
    {
        return false;
    }
}


void CushionOptimize::ObjectiveFunctionAndGradient(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
{
    TypeConvert type_converter;
    Eigen::VectorXd x_eigen;
    type_converter.AlglibArray2EigenVector(x, x_eigen);
    Eigen::VectorXd gradient;
    func = cross_section_info_->ObjectiveFunction2_WithGradient(x_eigen.data(), gradient);
    type_converter.EigenVector2AlglibArray(gradient, grad);
}
