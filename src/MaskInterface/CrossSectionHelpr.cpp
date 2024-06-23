#include "CrossSectionHelpr.h"

CrossSectionHelpr::CrossSectionHelpr()
{
}

CrossSectionHelpr::~CrossSectionHelpr()
{
}

double CrossSectionHelpr::TangentVector2Angle(const Eigen::Vector2d& _tangent) const
{
	// cos = tangent.x
	double angle = acos(_tangent.x());// in [0,pi]
	if (_tangent.y() < 0)
	{
		angle = -angle;
	}
	return angle;
}

void CrossSectionHelpr::TangentAngle2Vector(double _angle, Eigen::Vector2d& _tangent) const
{
	_tangent.x() = cos(_angle);
	_tangent.y() = sin(_angle);
}
