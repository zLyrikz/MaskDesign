#include "ConstraintFunctor.h"

ConstraintFunctor::ConstraintFunctor()
{
	cross_section_info_ = nullptr;
}

ConstraintFunctor::ConstraintFunctor(CrossSectionOptimizeInfo* _cross_section_info)
{
	cross_section_info_ = _cross_section_info;
}

ConstraintFunctor::~ConstraintFunctor()
{
}

bool ConstraintFunctor::operator()(const Eigen::VectorXd& _x) const
{
	return cross_section_info_->IsFeasible2(_x);
}
