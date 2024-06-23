#pragma once
#include<functional>

// https://en.wikipedia.org/wiki/Gaussian_quadrature
// with n points, exact result for polynomials of degree 2n-1 or less
// points and weights: https://pomax.github.io/bezierinfo/legendre-gauss.html
class GaussianQuadrature
{
public:
	GaussianQuadrature();
	~GaussianQuadrature();

	double TwoPointQuadrature(double _domain_start, double _domain_end, const std::function<double(double)>& _function) const;
};

