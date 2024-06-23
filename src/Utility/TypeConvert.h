#pragma once
#include "../Alglib/ap.h"
#include <Eigen/Core>

class TypeConvert
{
public:
	TypeConvert();
	~TypeConvert();
	// WARNING, when convert other type to Alglib, preset Alglib array lenth before use

	void AlglibArray2EigenVector(const alglib::real_1d_array& _x, Eigen::VectorXd& _eigen_x) const;
	void EigenVector2AlglibArray(const Eigen::VectorXd& _x, alglib::real_1d_array& _alglib_x) const;

	// WARNING allocate _x_data memory before use
	void AlglibArray2DataArray(const alglib::real_1d_array& _x, double* _x_data) const;
	void DataArray2AlglibArray(const double* _x_data, alglib::real_1d_array& _x) const;

};

