/*	This file LinearAlgebraInterface.h is part of scallop.
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

#ifndef SCALLOP_AUXILLARY_LINEARALGEBRAINTERFACE_H_
#define SCALLOP_AUXILLARY_LINEARALGEBRAINTERFACE_H_

#include <complex>
#include <vector>
#ifdef __MKL

#define MKL_INT int
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include <mkl_types.h>
#include <mkl_lapacke.h>
#include <mkl.h>

#else
extern "C"
{
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <cblas.h>
#include <lapacke.h>
}
#endif


namespace scallop
{
namespace auxillary
{

template<typename T>
class LinearAlgebraInterface
{
public:

	void matrix_times_diagonal_matrix(
			T const * matrix, size_t dim,
			T const * diagonalMatrix,
			T * resultMatrix) const;

	void matrix_times_diagonal_matrix(
			T * matrix, size_t dim,
			T const * diagonalMatrix) const;

	void matrix_times_matrix(
			T const * mleft, size_t dim,
			T const * mright,
			T * result) const;

	void call_gemm( bool transA, bool transB,
            int m, int n, int k,T alpha,T const * A, int lda,
            T const * B, int ldb, T beta,T* C,int ldc) const;

};

} /* namespace auxillary */
} /* namespace scallop */


#include "scallop/auxillary/src/LinearAlgebraInterface.hpp"
#endif /* SCALLOP_AUXILLARY_LINEARALGEBRAINTERFACE_H_ */
