/*	This file LinearAlgebraInterface.hpp is part of scallop.
 *
 *  scallop is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  scallop is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with scallop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Nov 3, 2016
 *      Author: A. Linscheid
 */

#include "scallop/auxillary/LinearAlgebraInterface.h"
#include <iostream>
#include <complex>

namespace scallop
{
namespace auxillary
{

template<typename T>
void LinearAlgebraInterface<T>::matrix_times_diagonal_matrix(
		typename std::vector<T>::const_iterator matrix, size_t dim,
		typename std::vector<T>::const_iterator diagonalMatrix,
		typename std::vector<T>::iterator resultMatrix) const
{
	auto itEnd = diagonalMatrix + dim;
	for ( ; diagonalMatrix != itEnd; ++diagonalMatrix )
		for ( size_t i=0 ; i<dim ; ++i )
		{
			*resultMatrix = (*diagonalMatrix)*(*matrix);
			++matrix;
			++resultMatrix;
		}
}

template<typename T>
void LinearAlgebraInterface<T>::matrix_times_matrix(
		typename std::vector<T>::const_iterator mleft, size_t dim,
		typename std::vector<T>::const_iterator mright,
		typename std::vector<T>::iterator result) const
{
	this->call_gemm(false,false,dim,dim,dim,T(1.0),&(*mleft),dim,&(*mright),dim,T(0.0),&(*result),dim);
}

} /* namespace auxillary */
} /* namespace scallop */
