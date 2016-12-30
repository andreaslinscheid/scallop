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

#include "scallop/auxillary/TypeMapComplex.h"
#include "scallop/auxillary/TemplateTypedefs.h"
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

	typedef typename auxillary::TypeMapComplex<T>::type bT;

	template<typename TD, typename TR>
	void matrix_times_diagonal_matrix(
			T const * matrix, size_t dim,
			TD const * diagonalMatrix,
			TR * resultMatrix) const;

	template<typename TD>
	void matrix_times_diagonal_matrix(
			T * matrix, size_t dim,
			TD const * diagonalMatrix) const;

	void matrix_times_matrix(
			T const * mleft, size_t dim,
			T const * mright,
			T * result) const;

	void call_gemm( bool transA, bool transB,
            int m, int n, int k,T alpha,T const * A, int lda,
            T const * B, int ldb, T beta,T* C,int ldc) const;

	void call_gemv(
			bool transA,
	        int m, int n, T alpha,T const * A, int lda,
	        T const * X, int incX, T beta,T * Y,int incY) const;

	/**
	 * 	Compute the eigensystem of the hermitian matrix.
	 *
	 *	We assume the data is located in upper(lower) triangle of the matrix if data_upper is true(false)
	 * 	Eigenvalues will be put into 'eigenval'
	 * 	NOTE: The matrix will be overwritten
	 */
	void hermitian_eigensystem(
			bool data_upper, bool comEV,
			T * matrix, int dim,
			bT * eigenval) const;

	void eigensystem(
			bool comEV,
			T * matrix, int dim,
			T * eigenval) const;

	/**
	 * Replace the square 'matrix' with its inverse.
	 *
	 * @param matrix					Row major ordered square matrix. D must implement data() and size() similar to vector.
	 * @param conditionNumber			If \p computeConditionNumber is true we return the condition number
	 * @param computeConditionNumber	Decide whether or not to compute condition number (by default we do)
	 */
	template<class D>
	void invert_square_matrix(D & matrix,
			bT & conditionNumber,
			bool computeConditionNumber = true) const;

	/**
	 * Replace the hermitian positive definite 'matrix' with its inverse.
	 *
	 * @param matrix					Row major ordered square matrix. D must implement data() and size() similar to vector.
	 */
	template<class D>
	void invert_positive_definite_hermitian_matrix(D & matrix) const;

	/**
	 * Check if the hermitian matrix 'matrix' is positive definite and if so replace it with its inverse or compute the eigenspectrum.
	 *
	 * @param inMatrix				Row major ordered input square matrix. D must implement data() and size() similar to vector.
	 * @param outmatrix				Row major ordered output square matrix. D must implement data() and size() similar to vector.
	 * 								Either the inverse matrix or the eigenvectors.
	 * @param is_pos_dev			set to true if the matrix turns out to be positive definite.
	 * @param eigenvalues			Significant if \p is_pos_dev is false. Then it contains the eigenvectors.
	 */
	template<class D, class V>
	void check_definite_invert_hermitian_matrix(D const& inMatrix,D & outMatrix, bool & is_pos_dev, V & eigenvalues) const;

private:

	mutable std::vector<int> IPIV;

	mutable typename auxillary::TemplateTypedefs<T>::scallop_vector workbuffer_;

	mutable typename auxillary::TemplateTypedefs<bT>::scallop_vector rWork_;

	template<class D>
	void determine_square_matrix_dim( D const& matrix, int & dim) const;

	void inversion_of_matrices_info_checks(int info, std::string what = "" ) const;

	int call_getri(int matrix_order, int n, T * a, int lda, const int * ipiv, T * work, int lwork) const;

	bT call_lange(int matrix_order, char norm, int m, int n, const T * a,  int lda) const;

	int call_getrf( int matrix_order, int m, int n, T * a, int lda, int * ipiv ) const;

	int call_gecon( int matrix_order, char norm, int n, const T* a, int lda, bT anorm, bT* rcond ) const;

	int call_heev( int matrix_order, char jobz, char uplo,
			   int n, T * a,
			   int lda, bT * w,
			   T * work, int lwork,
			   bT * rwork) const;

	int call_geev(  int matrix_order, char jobvl, char jobvr,
            int n, T* a,
            int lda, T* w,
            T* vl, int ldvl,
            T* vr, int ldvr,
            T* work, int lwork,
            bT* rwork) const;

	int call_potrf(
			int matrix_order, char uplo,int n,T * a,int lda ) const;

	int call_potri(
			int matrix_order, char uplo,int n,T * a,int lda ) const;
};

} /* namespace auxillary */
} /* namespace scallop */


#include "scallop/auxillary/src/LinearAlgebraInterface.hpp"
#endif /* SCALLOP_AUXILLARY_LINEARALGEBRAINTERFACE_H_ */
