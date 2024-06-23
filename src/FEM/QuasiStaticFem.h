#pragma once
#include "../MeshProcessing/TriangleMeshDistance.h"
#include "../Alglib/ap.h"
#include "../Alglib/optimization.h"
#include <functional>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Discregrid/All>


// for one 3D elasitc object 
// minimize elasitic energy; external force potential
// handle collision by SDF penalty of other objects
// not consider firction and self collision now
// in short: minimize objective function = E + W + Collision = total potential
// the unknown variable is x = displacement; so the displacement considered here is without those assigned in the boundary condition
class QuasiStaticFem
{
public:
	QuasiStaticFem();
	~QuasiStaticFem();

	void SetInitialConfiguration(const Eigen::VectorXd& _initial_nodes);
	void SetInitialDisplacement(const Eigen::VectorXd& _initial_displacement);
	void SetCollisionSdf(const tmd::TriangleMeshDistance* _other_objects, int _sdf_function_index = 0, const Discregrid::CubicLagrangeDiscreteGrid* _real_sdf = nullptr);
	void SetCollisionWeight(double _collision_weight);
	
	// set energy functions before optimization
	// Conjugate Gradient optimization
	void CgOptimize(alglib::ae_int_t _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);
	// L-BFGS optimization
	void LbfgsOptimize(alglib::ae_int_t _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1, bool _print_result = true);
	// L-BFGS penalty method (penalty term increases)
	void LbfgsPenaltyOptimize(alglib::ae_int_t _max_outer_iteration, double _max_collision, double _penalty_increase, 
		alglib::ae_int_t _max_inner_iteration, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1, int _M = 5,
		std::vector<Eigen::VectorXd>& _final_displacement_each_step = std::vector<Eigen::VectorXd>());

	// get the function pointer for setting those functions
	//std::function<double(const Eigen::VectorXd&)>* GetElasticEnergyFunction(); // for setting a energy function  
	//std::function<Eigen::VectorXd(const Eigen::VectorXd&)>* GetElasticEnergyGradient(); // for setting a energy gradiet function  
	std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)>* GetElasticEnergyValueAndGradient(); // for setting a energy term
	std::function<void(const Eigen::VectorXd&)>* GetInfoFunction();
	void GetNodeFromDisplacement(const Eigen::VectorXd& _displacement, Eigen::VectorXd& _nodes) const;
	void GetFinalDisplacement(Eigen::VectorXd& _final_displacement) const;

	// derivative info
	void FinalCollisionHessian(Eigen::SparseMatrix<double>& _hessian) const;

	// compute value and gradient in the same function; use common variables to save computation
	void CollisionValueAndGradient(const Eigen::VectorXd& _displacement, double& _value, Eigen::VectorXd& _gradient) const;

//public:
//	// test efficiency
//	double time_collision_ = 0.0;

private:
	void LbfgsOptimizeSetUp(alglib::minlbfgsstate& _state, int _M = 5, double _scale = 5, alglib::ae_int_t _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);


	// the objcetive function(total potential) in alglib data structure
	// this provide also the gradient info
	void ObjcetiveFunction_alglib(const alglib::real_1d_array& _displacement, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr) const;
	void ConstraintOptimizationFunctions_alglib(const alglib::real_1d_array& _displacement, alglib::real_1d_array& _functions, alglib::real_2d_array& _jacobian, void* ptr);

	// the gradient of the objective function(total potential)
	//void Gradient(const Eigen::VectorXd& _displacement, Eigen::VectorXd& _gradient) const;

	
	// sum ( min{sdf(x), 0} )^2
	// note: result adding the weight yet
	double CollisionPenalty(const Eigen::VectorXd& _displacement) const;// w.r.t. displacement
	//void CollisionGradient(const Eigen::VectorXd& _displacement, Eigen::VectorXd& _gradient) const;// gradient of the collision penalty; compute by finite difference

	double MinSdf0(const Eigen::Vector3d& _node) const;// min{sdf(x), 0}
	double Sdf(const Eigen::Vector3d& _node) const;// sdf(x)
	Eigen::Vector3d SdfGradient(const Eigen::Vector3d& _node) const;// D sdf(x)/Dx

	// Collision energy hessian w.r.t. unknown nodes
	// this is a 3x3 block diagonal matrix
	// not multiplied with weight
	void CollisionHessian(const Eigen::VectorXd& _displacement, Eigen::SparseMatrix<double>& _hessian) const;

private:
	// info of the problem
	int dof_ = 0;// the paramter degree of freedom(dimension) of the FEM problem
	int num_node_ = 0; // num_node_ = dof_ / 3 (since we consider objects in 3d space only)
	const int space_dim_ = 3;// = 3

	Eigen::VectorXd initial_nodes_;// initial configuration of the object; only the unknown ones
	Eigen::VectorXd initial_displacement_;// start point of the simulation
	Eigen::VectorXd final_displacement_;// final point of the simulation

	// elastic energy function by a displacment parameter (exclude the known, i.e. those in displacement boundary condition)
	//std::function<double(const Eigen::VectorXd&)> elastic_energy_;
	//std::function<Eigen::VectorXd(const Eigen::VectorXd&)> elastic_energy_gradient_;
	std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)> elastic_energy_value_and_gradient_;

	const tmd::TriangleMeshDistance* sdf_other_objects_ = nullptr;
	double collision_weight_ = 0.0;
	int sdf_function_index_ = 0;//for CubicLagrangeDiscreteGrid function query
	const Discregrid::CubicLagrangeDiscreteGrid* real_sdf_ = nullptr;
	std::function<double(const Eigen::VectorXd&)> collision_penalty_function_;//i.e. double CollisionPenalty(const Eigen::VectorXd& _displacement) const;
	std::function<double(const Eigen::Vector3d&)> min_sdf_0_function_;// i.e. double MinSdf0(const Eigen::Vector3d& _node)
	std::function<double(const Eigen::Vector3d&)> sdf_function_;//i.e. double Sdf(const Eigen::Vector3d& _node) const;
	std::function<Eigen::Vector3d(const Eigen::Vector3d&)> sdf_gradient_;

	// info during solution
	std::function<void(const Eigen::VectorXd&)> info_function_;// print element volume and energy information

};

