#pragma once
#include "CrossSectionOptimizeInfo.h"

class ObjectiveFunctor
{
public:
	ObjectiveFunctor();
	ObjectiveFunctor(CrossSectionOptimizeInfo* _cross_section_info);
	~ObjectiveFunctor();

	double operator()(const Eigen::VectorXd& _x, Eigen::VectorXd& _gradient, double _inexact_gradient_epsilon = 0.0) const;
private:

	CrossSectionOptimizeInfo* cross_section_info_;
};

