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

	fftw_complex * dataPtr = reinterpret_cast<fftw_complex*>(&data_[0]);

	int dim = static_cast<int>( spaceGrid_.size() );
	const parallel::MPIModule &mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_nproc() > 1 )
		dim--;

	int * inembed =NULL;
	int * onembed = NULL;
	int * n = new int [dim];
	int idist = 1;
	int odist = 1;
	int howmany = static_cast<int>( blockSize_*dimTimeFT_ );

	//Plan the R --> k transform
	if ( mpi.get_nproc() > 1 )
	{
		//plan Rx
		n[0] = spaceGrid_[0];

		int istride = static_cast<int>(blockSize_*dimTimeFT_*dimKProc_);
		int ostride = static_cast<int>(blockSize_*dimTimeFT_*dimKProc_);

		fftw3PlanKSingleFwd_ = fftw_plan_many_dft(
				1,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				-1,FFTW_PATIENT);

		// and plan the Ry and higher dimension
		for ( int i = 1; i < dim+1; i++)
			n[i] = spaceGrid_[i];

		istride = static_cast<int>(blockSize_*dimTimeFT_);
		ostride = static_cast<int>(blockSize_*dimTimeFT_);

		fftw3PlanKParaFwd_ = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				-1,FFTW_PATIENT);
	}
	else
	{
		for ( int i = 0; i < dim; i++)
			n[i] = spaceGrid_[i];

		int istride = static_cast<int>(blockSize_*dimTimeFT_);;
		int ostride = static_cast<int>(blockSize_*dimTimeFT_);;

		fftw3PlanKParaFwd_ = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				-1,FFTW_PATIENT);
	}

	//Plan the k --> R transform
	if ( mpi.get_nproc() > 1 )
	{
		//skip kx
		for ( int i = 1; i < dim+1; i++)
			n[i] = spaceGrid_[i];

		int istride = static_cast<int>(blockSize_*dimTimeFT_);
		int ostride = static_cast<int>(blockSize_*dimTimeFT_);

		fftw3PlanKParaBkwd_ = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				1,FFTW_PATIENT);

		//Now plan kx as the fastest running dimension
		n[0] = spaceGrid_[0];

		istride = static_cast<int>(blockSize_*dimTimeFT_*dimKProc_);
		ostride = static_cast<int>(blockSize_*dimTimeFT_*dimKProc_);

		fftw3PlanKSingleBkwd_ = fftw_plan_many_dft(
				1,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				1,FFTW_PATIENT);

		std::copy( dataCpy.begin(), dataCpy.end(), data_.begin() );
	}
	else
	{
		for ( int i = 0; i < dim; i++)
			n[i] = spaceGrid_[i];

		int istride = static_cast<int>(blockSize_*dimTimeFT_);;
		int ostride = static_cast<int>(blockSize_*dimTimeFT_);;

		fftw3PlanKParaBkwd_ = fftw_plan_many_dft(
				dim,n,howmany,
				dataPtr, inembed,istride, idist,
				dataPtr, onembed,ostride, odist,
				1,FFTW_PATIENT);
	}

	std::copy( dataCpy.begin(), dataCpy.end(), data_.begin() );

	delete [] n;
	delete [] inembed;
	delete [] onembed;
}

} /* namespace gw_flex */
} /* namespace scallop */
