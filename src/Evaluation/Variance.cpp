#include "Variance.h"
#include <numeric>
#include <iostream>
using std::cout;
using std::endl;
Variance::Variance()
{
}

Variance::~Variance()
{
}

double Variance::computeSum(const std::vector<double>& _values) const
{
	return std::accumulate(_values.begin(), _values.end(), 0.0);
}

double Variance::computeAverage(const std::vector<double>& _values) const
{
	double average = computeSum(_values) / double(_values.size());
	return average;
}

double Variance::computeVariance(const std::vector<double>& _values) const
{
	double average = computeAverage(_values);
	return computeVariance(average, _values);
}

double Variance::computeVariance(const std::vector<double>& _values, double& _sum) const
{
	_sum = computeSum(_values);
	double average = _sum / double(_values.size());

	return computeVariance(average, _values);
}

double Variance::computeVariance_WithGradient(const std::vector<double>& _values, Eigen::VectorXd& _gradient) const
{
	double sum = computeSum(_values);
	double average = sum / double(_values.size());
	double variance = computeVariance(average, _values);

	// gradient
	_gradient.setZero(_values.size());
	// term1: -2/n* (sum (x - average x)) (common term for all partial gradients)
	double term1 = 0.0;
	for (const double& iValue : _values)
	{
		term1 += iValue - average;
	}
	term1 *= -2.0 / double(_values.size());
	for (int iValue = 0; iValue < _values.size(); ++iValue)
	{
		_gradient(iValue) = term1 + 2.0 / double(_values.size()) * (_values[iValue] - average);
	}

	return variance;
}

double Variance::computeVariance(double _average, const std::vector<double>& _values) const
{
	double variance = 0.0;
	for (auto& i : _values)
	{
		double temp = (i - _average);
		variance += temp * temp;
	}

	variance /= double(_values.size());

	return variance;
}
