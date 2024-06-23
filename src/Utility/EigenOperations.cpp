#include "EigenOperations.h"

EigenOperations::EigenOperations()
{
}

EigenOperations::~EigenOperations()
{
}

void EigenOperations::ExtractSubmatrix_InputColumnMajor(const Eigen::SparseMatrix<double>& _matrix, 
	Eigen::SparseMatrix<double>& _submatrix, const std::vector<int>& _row_indices, const std::vector<int>& _col_indices)
{
	std::vector<Eigen::Triplet<double>> tripletList;
	for (int iCol = 0; iCol < _col_indices.size(); ++iCol)
	{
		if (_row_indices.size() > 0)
		{
			// at each col, iterate row of the submatrix from 0
			int iRow = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(_matrix, _col_indices[iCol]); it; ++it)
			{
				// at _col_indices[iCol]
				// for the current iRow, we want to find the matched nonzero row in matrix 
				// if (it.row() < _row_indices[iRow]) 
				//		++it until it.row()>=_row_indices[iRow]
				// 
				// now that it.row()>=_row_indices[iRow]
				//		++iRow until it.row()<=_row_indices[iRow]
				//		if ==
				//			we find you
				//		else <
				//			++it.row()

				if (it.row() >= _row_indices[iRow])
				{
					// find the row index that is not less than current iterate row: it.row()
					while (iRow < _row_indices.size() && _row_indices[iRow] < it.row())
					{
						iRow++;
					}
					if (iRow == _row_indices.size())
					{
						break;// end the current column search, move to next column
					}

					if (_row_indices[iRow] == it.row())
					{
						tripletList.push_back(Eigen::Triplet<double>(iRow, iCol, it.value()));

					}
				}
			}
		}
	}

	_submatrix.resize(_row_indices.size(), _col_indices.size());
	_submatrix.setFromTriplets(tripletList.begin(), tripletList.end());

}
