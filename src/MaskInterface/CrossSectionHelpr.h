#pragma once
#include <Eigen/Core>

// it's actually a helper
// common functions used in CrossSectionSurface and CrossSection
class CrossSectionHelpr
{
public:
	CrossSectionHelpr();
	~CrossSectionHelpr();

	// input tangent must be normalized
	double TangentVector2Angle(const Eigen::Vector2d& _tangent) const;
	void TangentAngle2Vector(double _angle, Eigen::Vector2d& _tangent) const;
};

