#pragma once
#include "problem.h"
#include "CrossSectionOptimizeInfo.h"

class ObjectiveForCppoptlib :
    public cppoptlib::Problem<double>
{
public:
    void setCrossSectionInfo(CrossSectionOptimizeInfo* cross_section_info);
    double value(const TVector& x);
    void gradient(const TVector& x, TVector& grad);
    double valueAndGradient(const TVector& x, TVector& grad);
    void saveResult(std::string _filepath, std::string _filename) const;

private:
    CrossSectionOptimizeInfo* cross_section_info_;
};

