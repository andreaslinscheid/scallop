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
		T const * matrix, size_t dim,
		T const * diagonalMatrix,
		T * resultMatrix) const
{
	for ( size_t i = 0; i < dim; ++i )
		for ( size_t j=0 ; j<dim ; ++j )
			resultMatrix[i*dim+j] = matrix[i*dim+j] *diagonalMatrix[i];
}

template<typename T>
void LinearAlgebraInterface<T>::matrix_times_diagonal_matrix(
		T * matrix, size_t dim,
		T const * diagonalMatrix) const
{
	for ( size_t i = 0; i < dim; ++i )
		for ( size_t j=0 ; j<dim ; ++j )
			matrix[i*dim+j] *= diagonalMatrix[i];
}

template<typename T>
void LinearAlgebraInterface<T>::matrix_times_matrix(
		T const * mleft, size_t dim,
		T const * mright,
		T * result) const
{
	this->call_gemm(false,false,dim,dim,dim,T(1.0),mleft,dim,mright,dim,T(0.0),result,dim);
}

} /* namespace auxillary */
} /* namespace scallop */
