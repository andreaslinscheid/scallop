/*	This file FFTBase.cpp is part of scallop.
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
 *  Created on: Oct 31, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/FFTBase.h"
#include "scallop/error_handling/error_handling.h"
#include <string>
#include <complex>

namespace scallop
{
namespace gw_flex
{

template<>
void
FFTBase< std::complex<double> >::plan_time_fft()
{
	auto dataCpy = data_;
	int dim = 1;
	int * n = new int [dim];
	n[0] = dimTimeFT_;
	int * inembed =NULL;
	int * onembed = NULL;
	size_t dataBlockPerK = blockSize_*dimTimeFT_;

	//create time plans for each k point
	fftw3PlanTimeFwd_ = new fftw_plan [dimKProc_];
	fftw3PlanTimeBkwd_ = new fftw_plan [dimKProc_];
	for ( size_t ik = 0 ; ik < dimKProc_; ++ik)
	{
		fftw_complex * dataPtr = reinterpret_cast<fftw_complex*>(&data_[ik*dataBlockPerK]);

		int howmany = static_cast<int>( blockSize_ );
		int idist = 1;
		int odist = 1;
		int istride = blockSize_;
		int ostride = blockSize_;

		fftw3PlanTimeFwd_[ik] = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				-1,FFTW_PATIENT);

		if ( ! fftw3PlanTimeFwd_[ik] )
			scallop::error_handling::Error(std::string()
				+"Unable to plan DFT for imaginary time to frequency at k: " + std::to_string(ik) );

		std::copy( dataCpy.begin(), dataCpy.end(), data_.begin() );

		fftw3PlanTimeBkwd_[ik] = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				1,FFTW_PATIENT);

		if ( ! fftw3PlanTimeFwd_[ik] )
			scallop::error_handling::Error(std::string()
				+"Unable to plan DFT for frequency to imaginary time at k: " + std::to_string(ik) );

		std::copy( dataCpy.begin(), dataCpy.end(), data_.begin() );
	}

	delete [] n;
	delete [] inembed;
	delete [] onembed;
}


template<>
void
FFTBase< std::complex<double> >::plan_space_fft()
{
	auto dataCpy = data_;

//	fftw_complex * dataPtr = reinterpret_cast<fftw_complex*>(&data_[0]);
//
//	int dim = static_cast<int>( dimKProc_ );
//	int * n = new int [dim];
//	int howmany = static_cast<int>( howmanyK_ );
//	for ( size_t i = 0; i < static_cast<size_t>(dim); i++)
//		n[i] = kgrid_.get_nk_dim(i);
//	int idist = static_cast<int>(kgrid_.get_nkpts());
//	int odist = static_cast<int>(kgrid_.get_nkpts());
//	int * inembed =NULL;
//	int * onembed = NULL;
//	int istride = 1;
//	int ostride = 1;
//
//	fftw3PlanFwd_ = fftw_plan_many_dft(
//			dim,n,howmany,
//			dataPtr, inembed,istride, idist,
//			dataPtr, onembed,ostride, odist,
//			-1,FFTW_PATIENT);
//
//	std::copy( data.begin(), data.end(), data_.begin() );
//
//	fftw3PlanBkwd_ = fftw_plan_many_dft(
//			dim,n,howmany,
//			dataPtr, inembed,istride, idist,
//			dataPtr, onembed,ostride, odist,
//			1,FFTW_PATIENT);
//
//	std::copy( data.begin(), data.end(), data_.begin() );
}

} /* namespace gw_flex */
} /* namespace scallop */
