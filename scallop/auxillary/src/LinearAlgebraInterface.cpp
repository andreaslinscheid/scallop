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
	cblas_cgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans,
			transB ? CblasTrans : CblasNoTrans, m, n, k,
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

	cblas_zgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans,
			transB ? CblasTrans : CblasNoTrans, m, n, k,
			reinterpret_cast<const void*>(&alpha), reinterpret_cast<const void*>(A), lda,
			reinterpret_cast<const void*>(B), ldb,
			reinterpret_cast<const void*>(&beta), reinterpret_cast<void*>(C), ldc);
}

template<>
void LinearAlgebraInterface< std::complex<double> >::call_gemv(
		bool transA,
        int m, int n, std::complex<double> alpha,std::complex<double> const * A, int lda,
        std::complex<double> const * X, int incX, std::complex<double> beta,std::complex<double> * Y,int incY) const
{
#ifdef DEBUG_BUILD
	auto c_to_ptr = [] ( std::complex<double> const * ptr ) {
		return reinterpret_cast<std::uintptr_t>(reinterpret_cast<const void*>(ptr));
	};

	//Check if all the pointers are at least 16 bit aligned
	assert( (c_to_ptr(A) % 16 ) == 0 );
	assert( (c_to_ptr(X) % 16 ) == 0 );
	assert( (c_to_ptr(Y) % 16 ) == 0 );
#endif

	cblas_zgemv(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, m, n,
			reinterpret_cast<const void*>(&alpha), reinterpret_cast<const void*>(A), lda,
			reinterpret_cast<const void*>(X), incX,
			reinterpret_cast<const void*>(&beta), reinterpret_cast<void*>(Y), incY);
}

template<>
int LinearAlgebraInterface< std::complex<float> >::call_getri(
		int matrix_order, int n, std::complex<float> * a, int lda,
		const int * ipiv, std::complex<float> * work, int lwork) const
{
	return LAPACKE_cgetri_work(matrix_order,n,a,lda,ipiv,work,lwork);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_getri(
		int matrix_order, int n, std::complex<double> * a, int lda,
		const int * ipiv, std::complex<double> * work, int lwork) const
{
	return LAPACKE_zgetri_work(matrix_order,n,a,lda,ipiv,work,lwork);
}

template<>
float LinearAlgebraInterface< std::complex<float> >::call_lange(
		int matrix_order, char norm, int m, int n, const std::complex<float> * a,  int lda) const
{
	return LAPACKE_clange(matrix_order,norm,m,n,a,lda);
}

template<>
double LinearAlgebraInterface< std::complex<double> >::call_lange(
		int matrix_order, char norm, int m, int n, const std::complex<double> * a,  int lda) const
{
	return LAPACKE_zlange(matrix_order,norm,m,n,a,lda);
}

template<>
int LinearAlgebraInterface< std::complex<float> >::call_getrf(
		int matrix_order, int m, int n, std::complex<float> * a, int lda, int * ipiv ) const
{
	return LAPACKE_cgetrf_work(matrix_order,m,n,a,lda,ipiv);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_getrf(
		int matrix_order, int m, int n, std::complex<double> * a, int lda, int * ipiv ) const
{
	return LAPACKE_zgetrf_work(matrix_order,m,n,a,lda,ipiv);
}

template<>
int LinearAlgebraInterface< std::complex<float> >::call_gecon(
		int matrix_order, char norm, int n, const std::complex<float> * a, int lda, float anorm, float * rcond ) const
{
	return LAPACKE_cgecon(matrix_order,norm, n, a, lda, anorm, rcond);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_gecon(
		int matrix_order, char norm, int n, const std::complex<double> * a, int lda, double anorm, double * rcond ) const
{
	return LAPACKE_zgecon(matrix_order,norm, n, a, lda, anorm, rcond);
}

template<>
int LinearAlgebraInterface< std::complex<float> >::call_heev(
		int matrix_order, char jobz, char uplo,
		   int n, std::complex<float>  * a,
		   int lda, float* w,
		   std::complex<float>  * work, int lwork,
		   float * rwork) const
{
	return LAPACKE_cheev_work(matrix_order,jobz, uplo, n, a, lda,w,work,lwork,rwork);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_heev(
		int matrix_order, char jobz, char uplo,
		   int n, std::complex<double>  * a,
		   int lda, double* w,
		   std::complex<double>  * work, int lwork,
		   double * rwork) const
{
	return LAPACKE_zheev_work(matrix_order,jobz, uplo, n, a, lda,w,work,lwork,rwork);
}

template<>
int LinearAlgebraInterface< std::complex<float> >::call_geev(
		int matrix_order, char jobvl, char jobvr,
        int n, std::complex<float>* a,
        int lda, std::complex<float>* w,
        std::complex<float>* vl, int ldvl,
        std::complex<float>* vr, int ldvr,
        std::complex<float>* work, int lwork,
        float* rwork) const
{
	return LAPACKE_cgeev_work(matrix_order,jobvl,jobvr,
	        n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork);
}

template<>
int LinearAlgebraInterface< std::complex<double> >::call_geev(
		int matrix_order, char jobvl, char jobvr,
        int n, std::complex<double>* a,
        int lda, std::complex<double>* w,
        std::complex<double>* vl, int ldvl,
        std::complex<double>* vr, int ldvr,
        std::complex<double>* work, int lwork,
        double* rwork) const
{
	return LAPACKE_zgeev_work(matrix_order,jobvl,jobvr,
	        n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork);
}

template <>
int LinearAlgebraInterface< std::complex<float> >::call_potrf(
		int matrix_order, char uplo,int n,std::complex<float> * a,int lda) const
{
	return LAPACKE_cpotrf(matrix_order,uplo,n,a,lda);
}

template <>
int LinearAlgebraInterface< std::complex<double> >::call_potrf(
		int matrix_order, char uplo,int n,std::complex<double> * a,int lda) const
{
	return LAPACKE_zpotrf(matrix_order,uplo,n,a,lda);
}

template <>
int LinearAlgebraInterface< std::complex<float> >::call_potri(
		int matrix_order, char uplo,int n,std::complex<float> * a,int lda ) const
{
	return LAPACKE_cpotri(matrix_order,uplo,n,a,lda);
}

template <>
int LinearAlgebraInterface< std::complex<double> >::call_potri(
		int matrix_order, char uplo,int n,std::complex<double> * a,int lda ) const
{
	return LAPACKE_zpotri(matrix_order,uplo,n,a,lda);
}

} /* namespace auxillary */
} /* namespace scallop */
