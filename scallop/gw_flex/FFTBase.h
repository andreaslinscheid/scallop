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

#include "scallop/parallel/GridDistribution.h"
#include "scallop/auxillary/AlignmentAllocator.h"
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

	FFTBase(FFTBase<T> const& other);

	FFTBase<T> & operator= (FFTBase<T> const& rhs);

	FFTBase<T> & operator= (FFTBase<T> && other );

	FFTBase(
			typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
			bool dataIsInKSpace,
			bool dataIsInTimeSpace,
			std::vector<size_t> spaceGridTotal,
			size_t dimTimeFT,
			size_t blockSize);

	~FFTBase();

	void initialize(
			typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
			bool dataIsInKSpace,
			bool dataIsInTimeSpace,
			std::vector<size_t> spaceGridTotal,
			size_t dimTimeFT,
			size_t blockSize);
	/**
	 * Provide read access to the data stored for a grid and time point.
	 *
	 * @param iGrid	Consecutive index of the processor local k or R grid.
	 * @param it	Index of time/frequency
	 * @return		Const Ptr to the beginning of the data block at this grid/time points.
	 */
	T const * read_phs_grid_ptr_block(size_t iGrid, size_t it ) const;

	/**
	 * Provide modifying access to the data stored for a grid and time point.
	 *
	 * @param iGrid	Consecutive index of the processor local k or R grid.
	 * @param it	Index of time/frequency
	 * @return		Ptr to the beginning of the data block at this grid/time points.
	 */
	T * write_phs_grid_ptr_block(size_t iGrid, size_t it );

	/**
	 * Provide read access to the data stored for given grid and time point.
	 *
	 * @param iGrid	Consecutive index of the data grid.
	 * @param it	Index of time/frequency
	 * @return		Const Ptr to the beginning of the data block at this grid/time points.
	 */
	T const * read_data_ptr_block(size_t iGrid, size_t it ) const;

	/**
	 * Provide modifying access to the data stored for given grid and time point.
	 *
	 * @param iGrid	Consecutive index of the data grid.
	 * @param it	Index of time/frequency
	 * @return		Ptr to the beginning of the data block at this grid/time points.
	 */
	T * write_data_ptr_block(size_t iGrid, size_t it );

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
	parallel::GridDistribution<T> const& get_spaceGrid_proc() const;

	void perform_space_fft();

	/**
	 * Flag that this object needs to be reinitialized.
	 *
	 * Sets isInit_ to false;
	 */
	void set_uninitialized();

	bool is_init() const;

	bool is_in_k_space() const;

	bool is_in_time_space() const;

	void set_time_space(bool newState);
protected:

	void perform_time_fft();

private:

	///keep track if this object is currently in reciprocal space
	bool isKSpace_ = true;

	///keep track if this object is currently in time space
	bool isTimeSpace_ = true;

	///Set to true once initialized
	bool isInit_ = false;

	///Description of the vectors in the grid
	parallel::GridDistribution<T> spaceGrid_;

	///Matsubara frequency dimension
	size_t dimTimeFT_ = 0;

	///chunk of data in the last dimension, for example Nambu-spin^2 * orbitals^2
	size_t blockSize_ = 0;

	///Plan the discrete Fourier transform part of the imaginary time to frequency transform
	///The convention with signs implies w_n==>\tau is forward.
	fftw_plan fftw3PlanTimeFwd_ = NULL;

	///Plan one discrete Fourier transform part of the frequency to imaginary time transform
	///	for each k grid point. The convention with signs implies \tau==>w_n is backward.
	fftw_plan fftw3PlanTimeBkwd_ = NULL;

	///Plan the forward Fourier transform in the non-parallelized K-grid dimensions
	///	The convention with signs implies R==>k is forward.
	fftw_plan fftw3PlanGridParaFwd_ = NULL;

	///Plan the backward Fourier transform in the non-parallelized K-grid dimensions
	///	The convention with signs implies k==>R is backward.
	fftw_plan fftw3PlanGridParaBkwd_ = NULL;

	///Plan the forward Fourier transforms in the parallelized K-grid dimension
	///for each	point in the remaining grid. The convention with signs implies R==>k is forward.
	fftw_plan fftw3PlanGridSingleFwd_ = NULL;

	///Plan the backward Fourier transform in the parallelized K-grid dimension
	///for each	point in the remaining grid. The convention with signs implies k==>R is backward.
	fftw_plan fftw3PlanGridSingleBkwd_ = NULL;

	///Main data container. Large (typically).
	typename auxillary::TemplateTypedefs<T>::scallop_vector data_;

	///Fourier buffer. We copy data here before FFTing to avoid having to keep track of complicated layouts.
	fftw_complex * FFTBuffer_ = NULL;

	void plan_time_fft();

	void plan_space_fft();

	void determine_fftbuffer_dim();

	void perform_R_to_k_fft();

	void perform_k_to_R_fft();

	void perform_time_to_freq_fft_in_R();

	void perform_freq_to_time_fft_in_R();

	void perform_time_to_freq_fft_in_k();

	void perform_freq_to_time_fft_in_k();

	void copy_or_fill_buffer_time_FFT( size_t indexBeginDataArray, bool fillBuffer );

	void copy_or_fill_buffer_first_FFT(size_t ikx, size_t indexTime, bool fillBuffer );

	void copy_or_fill_buffer_second_FFT(size_t ikyz, size_t indexTime, bool fillBuffer );

	void deallocate_fftw3_plans();
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/FFTBase.hpp"
#endif /* SCALLOP_GW_FLEX_FFTBASE_H_ */
