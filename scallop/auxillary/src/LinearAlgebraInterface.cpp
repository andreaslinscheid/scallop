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
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
			reinterpret_cast<const void*>(&alpha), reinterpret_cast<const void*>(A), lda,
			reinterpret_cast<const void*>(B), ldb,
			reinterpret_cast<const void*>(&beta), reinterpret_cast<void*>(C), ldc);
}

} /* namespace auxillary */
} /* namespace scallop */
