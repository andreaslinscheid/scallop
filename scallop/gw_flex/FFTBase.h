/*	This file FFTBase.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_FFTBASE_H_
#define SCALLOP_GW_FLEX_FFTBASE_H_

#include <fftw3.h>
#include <vector>

namespace scallop
{
namespace gw_flex
{

/**
 * 	This module handles the data and the Fourier transform in time and space.
 *
 * 	We implement the interface to the fftw3 library. Also the basic data storage happens here
 * 	because part of the FFTW3 plans is a data pointer. We like to keep these two things together.
 */
template<typename T>
class FFTBase
{
public:

	FFTBase();

	FFTBase(std::vector<T> data,
			std::vector<size_t> spaceGrid,
			size_t dimTimeFT,
			size_t blockSize);

	~FFTBase();

	void plan_ffts(
			std::vector<T> data,
			std::vector<size_t> spaceGrid,
			size_t dimTimeFT,
			size_t blockSize);

	typename std::vector<T>::const_iterator begin() const;

	typename std::vector<T>::const_iterator end() const;

	typename std::vector<T>::iterator data_begin_modify();

	typename std::vector<T>::iterator data_end_modify();

	size_t get_num_grid() const;

	size_t get_num_time() const;

	size_t get_data_block_size() const;

	/**
	 * Obtain the space grid this processor handles currently.
	 *
	 * @return Vector of dimension of space. The entries specify the current grid elements in that direction
	 * 			handled by this processor.
	 * 			Example: Total space dim = { 10 , 10 }, 2 processors.
	 * 					If this object is in k space : returns {5 , 10}
	 * 					If this object is in R space : returns {10, 5 }
	 */
	std::vector<size_t> get_spaceGrid_proc() const;

	void perform_R_to_k_fft();

	void perform_k_to_R_fft();

	bool is_init() const;

protected:

	void perform_time_to_freq_fft();

	void perform_freq_to_time_fft();

private:

	///Set to true once initialized
	bool isInit_ = false;

	///dimension of k-space minus the parallel x dimension
	size_t dimKProc_ = 0;

	///dimension of R-space minus the parallel y dimension
	size_t dimRProc_ = 0;

	///Number of k vectors in the grid
	std::vector<size_t> spaceGrid_;

	///Matsubara frequency dimension
	size_t dimTimeFT_ = 0;

	///chunk of data in the last dimension, for example Nambu-spin^2 * orbitals^2
	size_t blockSize_ = 0;

	///Plan one discrete Fourier transform part of the imaginary time to frequency transform
	///	for each k point. The convention with signs implies w_n==>\tau is forward.
	fftw_plan * fftw3PlanTimeFwd_ = NULL;

	//Plan one discrete Fourier transform part of the frequency to imaginary time transform
	//	for each k point. The convention with signs implies \tau==>w_n is backward.
	fftw_plan * fftw3PlanTimeBkwd_ = NULL;

	//Plan the forward Fourier transform in the non-parallelized K-grid dimensions
	//	The convention with signs implies R==>k is forward.
	fftw_plan fftw3PlanKParaFwd_ = NULL;

	//Plan the backward Fourier transform in the non-parallelized K-grid dimensions
	//	The convention with signs implies k==>R is backward.
	fftw_plan fftw3PlanKParaBkwd_ = NULL;

	//Plan the forward Fourier transform in the parallelized K-grid dimension
	//	 The convention with signs implies R==>k is forward.
	fftw_plan fftw3PlanKSingleFwd_ = NULL;

	//Plan the backward Fourier transform in the parallelized K-grid dimension
	//	 The convention with signs implies k==>R is backward.
	fftw_plan fftw3PlanKSingleBkwd_ = NULL;

	std::vector<T> data_;

	void plan_time_fft();

	void plan_space_fft();
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/FFTBase.hpp"
#endif /* SCALLOP_GW_FLEX_FFTBASE_H_ */
