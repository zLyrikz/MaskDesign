#include "ObjectiveForCppoptlib.h"
void ObjectiveForCppoptlib::setCrossSectionInfo(CrossSectionOptimizeInfo* cross_section_info)
{
    cross_section_info_ = cross_section_info;
}
double ObjectiveForCppoptlib::value(const TVector& x)
{
    TVector grad;
    return cross_section_info_->ObjectiveFunction2_WithGradient(x.data(), grad);
}

void ObjectiveForCppoptlib::gradient(const TVector& x, TVector& grad)
{
    cross_section_info_->ObjectiveFunction2_WithGradient(x.data(), grad);
}

double ObjectiveForCppoptlib::valueAndGradient(const TVector& x, TVector& grad)
{
    return cross_section_info_->ObjectiveFunction2_WithGradient(x.data(), grad);
}

void ObjectiveForCppoptlib::saveResult(std::string _filepath, std::string _filename) const
{
    return cross_section_info_->SaveCushionSurface(_filepath, _filename);
}
