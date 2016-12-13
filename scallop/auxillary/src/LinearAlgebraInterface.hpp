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
#include "scallop/error_handling/Error.h"
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

template<typename T>
void LinearAlgebraInterface<T>::hermitian_eigensystem(
		bool data_upper, bool comEV,
		T * matrix, int dim,
		bT * eigenval) const
{
	char v = ( comEV ?  'v' : 'n' );

	char dataLoc = ( data_upper ? 'U' : 'L' );

	int rworksize = ( 3*dim-2 > 1 ? 3 * dim - 2 : 1 ) ;
	if ( rWork_.size() != static_cast<size_t>(rworksize) )
		rWork_ = decltype(rWork_)(rworksize);

	T work_query;
	this->call_heev(
			LAPACK_ROW_MAJOR,
			v, dataLoc,
			dim, matrix, dim,
			eigenval,&work_query,-1,rWork_.data());
	size_t optimal_size = static_cast<size_t>( work_query.real() ) ;
	if ( workbuffer_.size() != optimal_size )
		workbuffer_ = decltype(workbuffer_)( optimal_size );


	int lwork = static_cast<int>(optimal_size);
	int info = this->call_heev(
			LAPACK_ROW_MAJOR,
			v, dataLoc,
			dim, matrix, dim,
			eigenval,
			workbuffer_.data(),
			lwork,
			rWork_.data());

	inversion_of_matrices_info_checks(info);

	//Ensure phase convention on eigenvectors (real part of the last component is non-negative)
	for ( int ib = 0 ; ib < dim ; ++ib )
	{
		if ( matrix[(dim-1)*dim + ib].real() < 0 )
			for ( int ibp = 0 ; ibp < dim ; ++ibp )
			{
				matrix[ibp*dim + ib] *= T(-1.0);
			}
	}
}

template<typename T>
template<class D>
void LinearAlgebraInterface<T>::invert_square_matrix(
		D & matrix,
		bT & conditionNumber,
		bool computeConditionNumber ) const
{
	int dim;
	this->determine_square_matrix_dim(matrix,dim);

	int info;
	if ( IPIV.size() != static_cast<size_t>(dim) )
	{
		IPIV.assign( dim, 0);
		T * work_query = new T [1];
		info = this->call_getri(LAPACK_ROW_MAJOR,dim,NULL,dim,NULL,work_query,-1);
		inversion_of_matrices_info_checks(info);
		workbuffer_.assign( static_cast<int>( work_query[0].real() ) , T(0) );
		delete [] work_query;
	}

	T * dptr = matrix.data();
	int * ipptr = IPIV.data();
	T * work = workbuffer_.data();
	int lwork = static_cast<int>(workbuffer_.size());

	bT matrixNorm;
	if ( computeConditionNumber )
		//Compute the matrix norm
		matrixNorm = this->call_lange(LAPACK_ROW_MAJOR,'1',dim,dim,dptr,dim);

	// L U factorization
	info = this->call_getrf(LAPACK_ROW_MAJOR,dim,dim,dptr,dim,ipptr);
	inversion_of_matrices_info_checks(info);

	if ( computeConditionNumber )
	{
		//Compute the condition number
		info = this->call_gecon(LAPACK_ROW_MAJOR,'1', dim, dptr, dim, matrixNorm, &conditionNumber);
		inversion_of_matrices_info_checks(info);
	}

	//compute the inverse
	info = this->call_getri(LAPACK_ROW_MAJOR,dim,dptr,dim,ipptr,work,lwork);
	inversion_of_matrices_info_checks(info);

}

template <typename T>
template<class D>
void LinearAlgebraInterface<T>::determine_square_matrix_dim( D const& matrix, int & dim) const
{
	dim = static_cast<int>( std::floor( std::sqrt( double(matrix.size()) ) + 0.5 ) );
	if ( static_cast<size_t>(dim*dim) != matrix.size() )
		error_handling::Error("Problem determining the size of the matrix. Re-think this strategy, maybe ...\n",4);
}

template <typename T>
void LinearAlgebraInterface<T>::inversion_of_matrices_info_checks( int info ) const
{
	if ( info == 0 )
		return;

	if ( info < 0 )
		std::cout << "WARNING input value # " << -info << " had an illegal value !" << std::endl;

	if ( info > 0 )
		std::cout << "WARNING LinAlg algorithm failed to converge: " << info << std::endl;
}

} /* namespace auxillary */
} /* namespace scallop */
