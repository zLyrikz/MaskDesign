#include "GaussianQuadrature.h"

GaussianQuadrature::GaussianQuadrature()
{
}

GaussianQuadrature::~GaussianQuadrature()
{
}

double GaussianQuadrature::TwoPointQuadrature(double _domain_start, double _domain_end, const std::function<double(double)>& _function) const
{
	double interval_change = (_domain_end - _domain_start) * 0.5;
	double translate = (_domain_end + _domain_start) * 0.5;

	const double point1 = -0.5773502691896257;
	const double point2 = 0.5773502691896257;

	double point1_translate = interval_change * point1 + translate;
	double point2_translate = interval_change * point2 + translate;

	//weight = 1.0

	double value = interval_change * (_function(point1_translate) + _function(point2_translate));

	return value;
}
