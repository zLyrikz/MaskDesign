#include "ChangeFrame.h"
#include <Eigen/LU>

ChangeFrame::ChangeFrame(const Eigen::Matrix4d& _frame1, const Eigen::Matrix4d& _frame2,
	Eigen::Matrix4d& _transformation)
{
	_transformation = Eigen::Matrix4d::Identity();
	Eigen::Matrix3d basis2_inverse = _frame2.block(0, 0, 3, 3).inverse();
	_transformation.block(0, 0, 3, 3) = basis2_inverse * _frame1.block(0, 0, 3, 3);
	_transformation.block(0, 3, 3, 1) = basis2_inverse * (_frame1.block(0, 3, 3, 1) - _frame2.block(0, 3, 3, 1));
}
