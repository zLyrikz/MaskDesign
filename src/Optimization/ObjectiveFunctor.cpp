#include "ObjectiveFunctor.h"

ObjectiveFunctor::ObjectiveFunctor()
{
	cross_section_info_ = nullptr;
}

ObjectiveFunctor::ObjectiveFunctor(CrossSectionOptimizeInfo* _cross_section_info)
{
	cross_section_info_ = _cross_section_info;
}

ObjectiveFunctor::~ObjectiveFunctor()
{
}

double ObjectiveFunctor::operator()(const Eigen::VectorXd& _x, Eigen::VectorXd& _gradient, double _inexact_gradient_epsilon) const
{
	return 	cross_section_info_->ObjectiveFunction2_WithGradient(_x.data(), _gradient, _inexact_gradient_epsilon);
}
