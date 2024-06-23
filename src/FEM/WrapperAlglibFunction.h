#pragma once
#include "../Alglib/ap.h"
#include <functional>

// for pass a member function pointer to a normal function pointer 

class WrapperAlglibFunction
{
public:

	//1 for QuasiStaticFem::ObjcetiveFunction_alglib
	static void wrapper_ObjcetiveFunction_alglib(const alglib::real_1d_array& _displacement, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr)
	{
		broker_ObjcetiveFunction_alglib(_displacement, _function_value, _gradient_value, _ptr);
	}
	// If error from inline, try use C++ 17 standard
	inline static std::function<
		void(const alglib::real_1d_array& _displacement, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr)>
		broker_ObjcetiveFunction_alglib; //must be inline

	//2 for QuasiStaticFem::ConstraintOptimizationFunctions_alglib
	static void wrapper_ConstraintOptimizationFunctions_alglib(const alglib::real_1d_array& _displacement, alglib::real_1d_array& _functions, alglib::real_2d_array& _jacobian, void* _ptr)
	{
		broker_ConstraintOptimizationFunctions_alglib(_displacement, _functions, _jacobian, _ptr);
	}
	// If error from inline, try use C++ 17 standard
	inline static std::function<
		void(const alglib::real_1d_array& _displacement, alglib::real_1d_array& _functions, alglib::real_2d_array& _jacobian, void* _ptr)>
		broker_ConstraintOptimizationFunctions_alglib; //must be inline

	//3 for CushionOptimize
	static void wrapper_CushionOptimizeFunctions(const alglib::real_1d_array& _x, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr)
	{
		broker_CushionOptimizeFunctions(_x, _function_value, _gradient_value, _ptr);
	}
	inline static std::function<
		void(const alglib::real_1d_array& _x, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr)>
		broker_CushionOptimizeFunctions; //must be inline
};

