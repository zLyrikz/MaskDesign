#include "NelderMeadRailOptimize.h"
#include "NelderMead.h"
#include <functional>
#include <iostream>
using std::cout;
using std::endl;

NelderMeadRailOptimize::NelderMeadRailOptimize()
{
	rail_info_ = nullptr;
}

NelderMeadRailOptimize::~NelderMeadRailOptimize()
{
}

void NelderMeadRailOptimize::SetRailInfoAssistant(RailOptimizeInfo* _rail_info)
{
	rail_info_ = _rail_info;
}

void NelderMeadRailOptimize::Optimize(Rail& _optimal_rail,
	double* _optimal_f, double _min_f_variance, double _step, int _convergence_check, int _max_f_evaluation,
	int* _count_f_evaluation, int* _num_restart)
{
	if (rail_info_ != nullptr)
	{
		std::function<double(const double*)> objective_function = std::bind(&RailOptimizeInfo::ObjectiveFunction, rail_info_, std::placeholders::_1);
		Optimize(_optimal_rail, objective_function, _optimal_f, _min_f_variance, _step, _convergence_check, _max_f_evaluation, _count_f_evaluation, _num_restart);
	}
}

void NelderMeadRailOptimize::Optimize(Rail& _optimal_rail, const std::function<double(const double*)>& _objective_function, 
	double* _optimal_f, double _min_f_variance, double _step, int _convergence_check, int _max_f_evaluation, int* _count_f_evaluation, int* _num_restart)
{

	NelderMead nelder_mead;

	const int dim_x = rail_info_->GetDimX();
	double* initial_x = new double[dim_x];
	for (int iDim = 0; iDim < dim_x; ++iDim)
	{
		initial_x[iDim] = rail_info_->GetInitialX()[iDim];
	}
	double* optimal_x = new double[dim_x];
	double* step = new double[dim_x];
	for (int iDim = 0; iDim < dim_x; ++iDim)
	{
		step[iDim] = _step;
	}

	double initial_f = _objective_function(initial_x);
	//cout << "initial f=" << initial_f << endl;

	nelder_mead.Optimize(_objective_function, dim_x, initial_x, optimal_x,
		_optimal_f, _min_f_variance, step, _convergence_check, _max_f_evaluation, _count_f_evaluation, _num_restart);
	nelder_mead.PrintFaultInfo();

	rail_info_->GetRailFromX(_optimal_rail, optimal_x);



	delete[] initial_x;
	delete[] optimal_x;
	delete[] step;
	
}
