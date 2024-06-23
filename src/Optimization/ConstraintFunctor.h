#pragma once
#include "CrossSectionOptimizeInfo.h"

class ConstraintFunctor
{
public:
	ConstraintFunctor();
	ConstraintFunctor(CrossSectionOptimizeInfo* _cross_section_info);
	~ConstraintFunctor();

	bool operator()(const Eigen::VectorXd& _x) const;

private:
	CrossSectionOptimizeInfo* cross_section_info_;
};

