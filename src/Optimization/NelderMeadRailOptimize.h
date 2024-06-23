#pragma once
#include "RailOptimizeInfo.h"

class NelderMeadRailOptimize
{
public:
	NelderMeadRailOptimize();
	~NelderMeadRailOptimize();

	void SetRailInfoAssistant(RailOptimizeInfo* _rail_info);



//    Output, double _optimal_f, the minimum value of the function.
//
//    Input, double _min_f_variance, the terminating limit for the variance
//    of function values.
//
//    Input, double _step[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int _convergence_check, the convergence check is carried out 
//    every _convergence_check iterations.
//
//    Input, int _max_f_evaluation, the maximum number of function 
//    evaluations.
//
//    Output, int *_count_f_evaluation, the number of function evaluations 
//    used.
//
//    Output, int *_num_restart, the number of restarts.
//
//    Output, int *_ifault, error indicator.
//    0, no errors detected.
//    1, _min_f_variance, N, or _convergence_check has an illegal value.
//    2, iteration terminated because _max_f_evaluation was exceeded without convergence.

	void Optimize(Rail& _optimal_rail,
		double* _optimal_f, double _min_f_variance, double _step, int _convergence_check, int _max_f_evaluation,
		int* _count_f_evaluation, int* _num_restart);

private:
	void Optimize(Rail& _optimal_rail, const std::function<double(const double*)>& _objective_function,
		double* _optimal_f, double _min_f_variance, double _step, int _convergence_check, int _max_f_evaluation,
		int* _count_f_evaluation, int* _num_restart);
private:
	// for provide some rail information for optimization
	RailOptimizeInfo* rail_info_;
};

