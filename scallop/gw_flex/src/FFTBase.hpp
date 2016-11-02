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
		size_t dimKProc,
		std::vector<size_t> howmanyK,
		size_t dimTimeFT,
		size_t blockSize)
{
	this->plan_ffts(
			std::move(data),dimKProc,howmanyK,dimTimeFT,blockSize);
}

template<typename T>
void
FFTBase<T>::plan_ffts(
		std::vector<T> data,
		size_t dimKProc,
		std::vector<size_t> howmanyK,
		size_t dimTimeFT,
		size_t blockSize)
{
	dimKProc_ = dimKProc;
	howmanyK_ = howmanyK;
	dimTimeFT_ = dimTimeFT;
	blockSize_ = blockSize;
	data_ = std::move(data);

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

} /* namespace gw_flex */
} /* namespace scallop */
