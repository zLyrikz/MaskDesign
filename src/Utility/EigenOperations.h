#pragma once
#include <vector>
#include <Eigen/SparseCore>

class EigenOperations
{
public:
	EigenOperations();
	~EigenOperations();

	// require row indices in increasing order
	static void ExtractSubmatrix_InputColumnMajor(const Eigen::SparseMatrix<double>& _matrix,
		Eigen::SparseMatrix<double>& _submatrix, const std::vector<int>& _row_indices, const std::vector<int>& _col_indices);
};

