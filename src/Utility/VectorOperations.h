#pragma once
#include <vector>
class VectorOperations
{
public:
	VectorOperations();

	// reference about template and inline https://www.cnblogs.com/stemon/p/4593582.html
	template<class T>
	void AppendBack(std::vector<T>& _origin, const std::vector<T>& _addition);

	template<class T>
	void AppendFront(std::vector<T>& _origin, const std::vector<T>& _addition);

	// 3D point like Mesh::Point, Eigen::Vector3
	// unique vector by pushback duplicate elements to _unique only once
	// also return the choosen element id
	template<class Point3D>
	void UniquePoint(const std::vector<Point3D>& _old, std::vector<Point3D>& _unique, std::vector<int>& _element_id, double _absoulte_error = 1e-5);
};

template<class T>
inline void VectorOperations::AppendBack(std::vector<T>& _origin, const std::vector<T>& _addition)
{
	// some references
	// right value reference https://blog.csdn.net/xiaolewennofollow/article/details/52559306
	// std::move https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used

	_origin.reserve(_origin.size() + _addition.size());
	_origin.insert(_origin.end(), _addition.begin(), _addition.end());
}

template<class T>
inline void VectorOperations::AppendFront(std::vector<T>& _origin, const std::vector<T>& _addition)
{
	_origin.reserve(_origin.size() + _addition.size());
	_origin.insert(_origin.begin(), _addition.begin(), _addition.end());
}

template<class Point3D>
inline void VectorOperations::UniquePoint(const std::vector<Point3D>& _old, std::vector<Point3D>& _unique, std::vector<int>& _element_id, double _absoulte_error)
{
	if (_old.size() > 0)
	{
		_unique.push_back(_old[0]);
		_element_id.push_back(0);
		for (int iCurrent = 1; iCurrent < _old.size(); iCurrent++)
		{
			// check if no previous element == current elemnt
			bool tag = true;
			for (int iPrevious = 0; iPrevious < iCurrent; iPrevious++)
			{
				if (abs(_old[iPrevious][0] - _old[iCurrent][0]) < _absoulte_error &&
					abs(_old[iPrevious][1] - _old[iCurrent][1]) < _absoulte_error &&
					abs(_old[iPrevious][2] - _old[iCurrent][2]) < _absoulte_error)
				{
					tag = false;// does exist
					break;
				}
			}
			if (tag == true)
			{
				_unique.push_back(_old[iCurrent]);
				_element_id.push_back(iCurrent);
			}

		}
	}

}
