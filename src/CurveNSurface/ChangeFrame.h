#pragma once
#include <Eigen/Core>

/*
* find relation matrix of coordinate X1, X2 of a point under two frames: [O; b1, b2, b3], [OO; bb1, bb2, bb3]
* 
* frame:
* 	// (b1x b2x b3x | Ox)
	// (b1y b2y b3z | Oy)
	// (b1z b2z b3z | Oz)
	// -------------------
	// (0   0   0   |  1)
* 
* transformation:
* 	// (R R R | T)
	// (R R R | T)
	// (R R R | T)
	// ----------
	// (0 0 0 | 1)
*
* then, coordinate of second frame X2 = transforamtion * X1
*/

class ChangeFrame
{
public:
	// change frame from 1 to 2
	// _frame1 * X1 = _frame2 * X2
	// X2 = transforamtion * X1
	ChangeFrame(const Eigen::Matrix4d& _frame1, const Eigen::Matrix4d& _frame2,
		Eigen::Matrix4d& _transformation); 


};

