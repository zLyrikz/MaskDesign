#include "TypeConvert.h"
#include <iostream>
using std::cout;
using std::endl;
TypeConvert::TypeConvert()
{
}

TypeConvert::~TypeConvert()
{
}

void TypeConvert::AlglibArray2EigenVector(const alglib::real_1d_array& _x, Eigen::VectorXd& _eigen_x) const
{
	int size = _x.length();
	_eigen_x.resize(size);
	for (int i = 0; i < size; ++i)
	{
		_eigen_x[i] = _x[i];
	}
}

void TypeConvert::EigenVector2AlglibArray(const Eigen::VectorXd& _x, alglib::real_1d_array& _alglib_x) const
{
	int size = _x.rows();

	for (int i = 0; i < size; ++i)
	{
		_alglib_x[i] = _x[i];
	}

}

void TypeConvert::AlglibArray2DataArray(const alglib::real_1d_array& _x, double* _x_data) const
{
	int size = _x.length();
	for (int i = 0; i < size; ++i)
	{
		_x_data[i] = _x[i];
	}
}

void TypeConvert::DataArray2AlglibArray(const double* _x_data, alglib::real_1d_array& _x) const
{
	int size = _x.length();
	for (int i = 0; i < size; ++i)
	{
		_x[i] = _x_data[i];
	}
}
