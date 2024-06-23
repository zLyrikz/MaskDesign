#pragma once
#include <vector>
#include <Eigen\Core>

class Variance
{
public:
	Variance();
	~Variance();

public:
	double computeSum(const std::vector<double>& _values) const;
	double computeAverage(const std::vector<double>& _values) const;
	double computeVariance(const std::vector<double>& _values) const;
	double computeVariance(const std::vector<double>& _values, double& _sum) const;
	double computeVariance_WithGradient(const std::vector<double>& _values, Eigen::VectorXd& _gradient) const;
private:
	double computeVariance(double _average, const std::vector<double>& _values) const;
};

