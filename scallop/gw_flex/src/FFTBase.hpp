/*	This file FFTBase.hpp is part of scallop.
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
#include "scallop/parallel/MPIModule.h"
#include <fftw3.h>

namespace scallop
{
namespace gw_flex
{

template<typename T>
FFTBase<T>::~FFTBase()
{
	for ( size_t ik = 0 ; ik < dimKProc_ ; ++ik)
	{
		fftw_destroy_plan(fftw3PlanTimeFwd_[ik]);
		fftw_destroy_plan(fftw3PlanTimeBkwd_[ik]);
	}
	delete [] fftw3PlanTimeFwd_;
	delete [] fftw3PlanTimeBkwd_;
 	fftw_destroy_plan(fftw3PlanKParaFwd_);
	fftw_destroy_plan(fftw3PlanKParaBkwd_);
	fftw_destroy_plan(fftw3PlanKSingleFwd_);
	fftw_destroy_plan(fftw3PlanKSingleBkwd_);
}

template<typename T>
FFTBase<T>::FFTBase()
{

}

template<typename T>
FFTBase<T>::FFTBase(std::vector<T> data,
		std::vector<size_t> spaceGrid,
		size_t dimTimeFT,
		size_t blockSize)
{
	this->plan_ffts(
			std::move(data),spaceGrid,dimTimeFT,blockSize);
}

template<typename T>
void
FFTBase<T>::plan_ffts(
		std::vector<T> data,
		std::vector<size_t> spaceGrid,
		size_t dimTimeFT,
		size_t blockSize)
{
	isInit_ = true;
	spaceGrid_ = std::move(spaceGrid);
	dimTimeFT_ = dimTimeFT;
	blockSize_ = blockSize;
	data_ = std::move(data);

	const parallel::MPIModule &mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_nproc() > 1 )
	{
		if ( spaceGrid_.size() < 2 )
			error_handling::Error("No parallelization for D < 2");

		//We parallelize along the x direction in k space and the y direction in real space.
		dimKProc_ = 1;
		for ( size_t i = 1 ; i < spaceGrid_.size(); i++ )
			dimKProc_ *= spaceGrid_[i];

		dimRProc_ = 1;
		for ( size_t i = 0 ; i < spaceGrid_.size(); i++ )
			if ( i != 1 )
				dimRProc_ *= spaceGrid_[i];

		if ( spaceGrid_[0] !=  spaceGrid_[1] )
			error_handling::Error("Because of parallelization issues, "
					"no x /= y grid dimension is supported in parallel at present");
	}
	else
	{
		dimKProc_ = 1;
		for ( auto nki : spaceGrid )
			dimKProc_ *= nki;
		dimRProc_ = dimKProc_;
	}

	this->plan_time_fft();

	this->plan_space_fft();
}

template<typename T>
size_t FFTBase<T>::get_num_grid() const
{
	return dimKProc_;
}

template<typename T>
size_t FFTBase<T>::get_num_time() const
{
	return dimTimeFT_;
}

template<typename T>
size_t FFTBase<T>::get_data_block_size() const
{
	return blockSize_;
}

template<typename T>
std::vector<size_t> FFTBase<T>::get_spaceGrid_proc() const
{
	return spaceGrid_;
}

template<typename T>
typename std::vector<T>::iterator
FFTBase<T>::data_begin_modify()
{
	return data_.begin();
}

template<typename T>
typename std::vector<T>::iterator
FFTBase<T>::data_end_modify()
{
	return data_.end();
}

template<typename T>
typename std::vector<T>::const_iterator
FFTBase<T>::begin() const
{
	return data_.begin();
}

template<typename T>
typename std::vector<T>::const_iterator
FFTBase<T>::end() const
{
	return data_.end();
}

template<typename T>
void FFTBase<T>::perform_time_to_freq_fft()
{
	for ( size_t ik = 0 ; ik < dimKProc_; ++ik)
		fftw_execute( fftw3PlanTimeBkwd_[ik] );
}

template<typename T>
void FFTBase<T>::perform_freq_to_time_fft()
{
	for ( size_t ik = 0 ; ik < dimKProc_; ++ik)
		fftw_execute( fftw3PlanTimeFwd_[ik] );
}

template<typename T>
void FFTBase<T>::perform_R_to_k_fft()
{
	const parallel::MPIModule &mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_nproc() == 1 )
	{
		fftw_execute( fftw3PlanKParaFwd_ );
	}
	else // ( mpi.get_nproc() == 1 )
	{

	}

}

template<typename T>
void FFTBase<T>::perform_k_to_R_fft()
{
	size_t spaceGridTotal = 1;
	for ( auto kd : spaceGrid_ )
		spaceGridTotal *= kd;

	const parallel::MPIModule &mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_nproc() == 1 )
	{
		fftw_execute( fftw3PlanKParaBkwd_ );
	}
	else // ( mpi.get_nproc() == 1 )
	{

	}

	//Normalization
	for ( auto &d : data_ )
		d *= 1.0/spaceGridTotal;
}

template<typename T>
bool FFTBase<T>::is_init() const
{
	return isInit_;
}

} /* namespace gw_flex */
} /* namespace scallop */
