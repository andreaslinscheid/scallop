/*	This file LinearAlgebraInterface.cpp is part of scallop.
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
#include <assert.h>
#include <cstdint>

namespace scallop
{
namespace auxillary
{


template<>
void LinearAlgebraInterface< std::complex<float> >::call_gemm(
		bool transA, bool transB,
        int m, int n, int k,std::complex<float> alpha,std::complex<float> const * A, int lda,
        std::complex<float> const * B, int ldb, std::complex<float> beta,std::complex<float> * C,int ldc) const
{
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
			reinterpret_cast<const void*>(&alpha), reinterpret_cast<const void*>(A), lda,
			reinterpret_cast<const void*>(B), ldb,
			reinterpret_cast<const void*>(&beta), reinterpret_cast<void*>(C), ldc);
}

template<>
void LinearAlgebraInterface< std::complex<double> >::call_gemm(
		bool transA, bool transB,
        int m, int n, int k,std::complex<double> alpha,std::complex<double> const * A, int lda,
        std::complex<double> const * B, int ldb, std::complex<double> beta,std::complex<double> * C,int ldc) const
{
#ifdef DEBUG_BUILD
	auto c_to_ptr = [] ( std::complex<double> const * ptr ) {
		return reinterpret_cast<std::uintptr_t>(reinterpret_cast<const void*>(ptr));
	};

	//Check if all the pointers are at least 16 bit aligned
	assert( (c_to_ptr(A) % 16 ) == 0 );
	assert( (c_to_ptr(B) % 16 ) == 0 );
	assert( (c_to_ptr(C) % 16 ) == 0 );
#endif

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
			reinterpret_cast<const void*>(&alpha), reinterpret_cast<const void*>(A), lda,
			reinterpret_cast<const void*>(B), ldb,
			reinterpret_cast<const void*>(&beta), reinterpret_cast<void*>(C), ldc);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_getri(
		int matrix_order, int n, std::complex<double> * a, int lda,
		const int * ipiv, std::complex<double> * work, int lwork) const
{
	return LAPACKE_zgetri_work(matrix_order,n,a,lda,ipiv,work,lwork);
}

template<>
double LinearAlgebraInterface< std::complex<double> >::call_lange(
		int matrix_order, char norm, int m, int n, const std::complex<double> * a,  int lda) const
{
	return LAPACKE_zlange(matrix_order,norm,m,n,a,lda);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_getrf(
		int matrix_order, int m, int n, std::complex<double> * a, int lda, int * ipiv ) const
{
	return LAPACKE_zgetrf_work(matrix_order,m,n,a,lda,ipiv);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_gecon(
		int matrix_order, char norm, int n, const std::complex<double> * a, int lda, double anorm, double * rcond ) const
{
	return LAPACKE_zgecon(matrix_order,norm, n, a, lda, anorm, rcond);
}

} /* namespace auxillary */
} /* namespace scallop */
