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
#include <assert.h>

namespace scallop
{
namespace gw_flex
{

template<typename T>
FFTBase<T> & FFTBase<T>::operator= (FFTBase<T> && other )
{
	assert(this != &other);
	isKSpace_ 					= std::move(other.isKSpace_);
	isTimeSpace_ 				= std::move(other.isTimeSpace_);
	isInit_ 					= std::move(other.isInit_);
	spaceGrid_					= std::move(other.spaceGrid_);
	dimTimeFT_ 					= std::move(other.dimTimeFT_);
	blockSize_ 					= std::move(other.blockSize_);
	fftw3PlanTimeFwd_ 			= std::move(other.fftw3PlanTimeFwd_);
	other.fftw3PlanTimeFwd_ = NULL;
	fftw3PlanTimeBkwd_ 			= std::move(other.fftw3PlanTimeBkwd_);
	other.fftw3PlanTimeBkwd_ = NULL;
	fftw3PlanGridParaFwd_		= std::move(other.fftw3PlanGridParaFwd_);
	other.fftw3PlanGridParaFwd_ = NULL;
	fftw3PlanGridParaBkwd_ 		= std::move(other.fftw3PlanGridParaBkwd_);
	other.fftw3PlanGridParaBkwd_ = NULL;
	fftw3PlanGridSingleFwd_ 	= std::move(other.fftw3PlanGridSingleFwd_);
	other.fftw3PlanGridSingleFwd_ = NULL;
	fftw3PlanGridSingleBkwd_ 	= std::move(other.fftw3PlanGridSingleBkwd_);
	other.fftw3PlanGridSingleBkwd_ = NULL;
	data_						= std::move(other.data_);
	FFTBuffer_ 					= std::move(other.FFTBuffer_);
	other.FFTBuffer_ = NULL;
	return *this;
}

template<typename T>
FFTBase<T>::~FFTBase()
{
	this->deallocate_fftw3_plans();
}

template<typename T>
void FFTBase<T>::deallocate_fftw3_plans()
{
	//They will never be allocated alone
	if ( fftw3PlanTimeFwd_ )
	{
		fftw_destroy_plan(fftw3PlanTimeFwd_);
		fftw3PlanTimeFwd_ = NULL;
		fftw_destroy_plan(fftw3PlanTimeBkwd_);
		fftw3PlanTimeBkwd_ = NULL;
		fftw_destroy_plan(fftw3PlanGridParaFwd_);
		fftw3PlanGridParaFwd_ = NULL;
		fftw_destroy_plan(fftw3PlanGridParaBkwd_);
		fftw3PlanGridParaBkwd_ = NULL;
		fftw_destroy_plan(fftw3PlanGridSingleFwd_);
		fftw3PlanGridSingleFwd_ = NULL;
		fftw_destroy_plan(fftw3PlanGridSingleBkwd_);
		fftw3PlanGridSingleBkwd_ = NULL;
		fftw_free( FFTBuffer_ );
		FFTBuffer_ = NULL;
	}
}

template<typename T>
FFTBase<T>::FFTBase()
{

}

template<typename T>
FFTBase<T>::FFTBase(FFTBase<T> const& other)
{
	*this = other;
}

template<typename T>
FFTBase<T> & FFTBase<T>::operator= (FFTBase<T> const& rhs)
{
	if (this != &rhs)
	{
		this->deallocate_fftw3_plans();

		isKSpace_ 	= rhs.isKSpace_;
		isTimeSpace_= rhs.isTimeSpace_;
		isInit_ 	= rhs.isInit_;
		spaceGrid_ 	= rhs.spaceGrid_;
		dimTimeFT_ 	= rhs.dimTimeFT_;
		blockSize_ 	= rhs.blockSize_;
		data_ 		= rhs.data_;

		this->determine_fftbuffer_dim();

		this->plan_time_fft();

		this->plan_space_fft();
	}
	return *this;
}

template<typename T>
FFTBase<T>::FFTBase(
		typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
		bool dataIsInKSpace,
		bool dataIsInTimeSpace,
		std::vector<size_t> spaceGridTotal,
		size_t dimTimeFT,
		size_t blockSize)
{
	this->initialize(dataIsInKSpace,
			std::move(data),spaceGridTotal,dimTimeFT,blockSize);
}

template<typename T>
void
FFTBase<T>::initialize(
		typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
		bool dataIsInKSpace,
		bool dataIsInTimeSpace,
		std::vector<size_t> spaceGridTotal,
		size_t dimTimeFT,
		size_t blockSize)
{
	this->deallocate_fftw3_plans();

	isInit_ = true;
	isKSpace_ 	= dataIsInKSpace;
	isTimeSpace_= dataIsInTimeSpace;
	spaceGrid_.distribute_grid( std::move(spaceGridTotal) );
	dimTimeFT_ = dimTimeFT;
	blockSize_ = blockSize;

	//The data passed to this routine only has elements for
	//physical grid points, not for the zero padded ones.
	size_t numDataPts = spaceGrid_.get_num_grid_data()*dimTimeFT_*blockSize_;
	data_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( numDataPts  );

	if ( ! data.empty() )
	{
		size_t nG = dataIsInKSpace ? spaceGrid_.get_num_k_grid() : spaceGrid_.get_num_R_grid();
		for ( size_t ig = 0 ; ig < nG; ++ig )
		{
			size_t sizePerGridPoint = dimTimeFT_*blockSize_;
			size_t internalDataGridIndex =  dataIsInKSpace ?
					spaceGrid_.k_conseq_local_to_data_conseq( ig )
					: spaceGrid_.R_conseq_local_to_data_conseq( ig );

			for ( size_t i = 0 ; i < sizePerGridPoint; ++i)
				data_[internalDataGridIndex*sizePerGridPoint+i] = data[ig*sizePerGridPoint+i];
		}
	}

	this->determine_fftbuffer_dim();

	this->plan_time_fft();

	this->plan_space_fft();
}

template<typename T>
void FFTBase<T>::determine_fftbuffer_dim()
{
	auto max = [] (size_t i, size_t j) { return i < j ? j :i ; };

	//Check the time dim, we do all FFT in time for one grid point at once
	size_t fftBufferDim = dimTimeFT_*blockSize_;

	//Check the first k to R space transform; we do all FFT for the blockSize_ at once
	//Note that k is parallelized in the kx direction. Thus our first FFT is the y (and possibly z)
	size_t nEleKToRFirst = spaceGrid_.get_num_k_grid()/spaceGrid_.get_k_grid().front();
	fftBufferDim = max(fftBufferDim,nEleKToRFirst*blockSize_);

	//Check the second k to R space transform; we do all FFT for the blockSize_ at once
	//Our second FFT is the kx direction
	size_t nEleKToRSecond = spaceGrid_.get_grid().front();
	fftBufferDim = max(fftBufferDim,nEleKToRSecond*blockSize_);

	//Check the first R to k space transform; we do all FFT for the blockSize_ at once
	//Note that R is parallelized in the last direction.
	size_t nEleRTokFirst = spaceGrid_.get_grid().front();
	fftBufferDim = max(fftBufferDim,nEleRTokFirst*blockSize_);

	//Check the second R to k space transform; we do all FFT for the blockSize_ at once
	//Note that R is parallelized in the kz direction. Thus our FFT is the x (and possibly y)
	size_t nEleRTokSecond = spaceGrid_.get_num_R_grid()/spaceGrid_.get_R_grid().back();
	fftBufferDim = max(fftBufferDim,nEleRTokSecond*blockSize_);

	FFTBuffer_ = fftw_alloc_complex( fftBufferDim );
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
parallel::GridDistribution<T> const& FFTBase<T>::get_spaceGrid_proc() const
{
	return spaceGrid_;
}

template<typename T>
void FFTBase<T>::copy_or_fill_buffer_time_FFT( size_t indexBeginDataArray, bool fillBuffer )
{
	auto time_FFT_layout = [&] ( size_t itime, size_t ib )
		{
			T * bufferAsCmplxArray = reinterpret_cast<T*>(FFTBuffer_);
			T * thisLocation = bufferAsCmplxArray + (ib*dimTimeFT_ + itime);
			return thisLocation;
		};

	auto it = data_.begin();
	if ( fillBuffer )
	{
		for ( size_t itime = 0 ; itime < dimTimeFT_ ; ++itime )
			for ( size_t ib = 0 ; ib < blockSize_ ; ++ib )
			{
				size_t index = indexBeginDataArray + itime*blockSize_+ib;
				( *time_FFT_layout(itime,ib) ) = *(it+index);
			}
	}
	else
	{
		for ( size_t itime = 0 ; itime < dimTimeFT_ ; ++itime )
			for ( size_t ib = 0 ; ib < blockSize_ ; ++ib )
			{
				size_t index = indexBeginDataArray + itime*blockSize_+ib;
				( *(it+index) ) = ( *time_FFT_layout(itime,ib) );
			}
	}
}

template<typename T>
void FFTBase<T>::perform_time_to_freq_fft_in_k()
{
	for ( size_t ik = 0 ; ik < spaceGrid_.get_num_k_grid(); ++ik)
	{
		//copy data to the buffer such that for each block,
		//the time dimension is a linear array
		size_t dataGridIndex = spaceGrid_.k_conseq_local_to_data_conseq( ik );
		size_t thisFFTDataBegin = dataGridIndex*blockSize_*dimTimeFT_;
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,true);

		//FFT
		fftw_execute( fftw3PlanTimeBkwd_ );

		//copy back
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,false);
	}
}

template<typename T>
void FFTBase<T>::perform_freq_to_time_fft_in_k()
{
	for ( size_t ik = 0 ; ik < spaceGrid_.get_num_k_grid(); ++ik)
	{
		//copy data to the buffer such that for each block,
		//the time dimension is a linear array
		size_t dataGridIndex = spaceGrid_.k_conseq_local_to_data_conseq( ik );
		size_t thisFFTDataBegin = dataGridIndex*blockSize_*dimTimeFT_;
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,true);

		//FFT
		fftw_execute( fftw3PlanTimeFwd_ );

		//copy back
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,false);
	}
}

template<typename T>
void FFTBase<T>::perform_time_to_freq_fft_in_R()
{
	for ( size_t iR = 0 ; iR < spaceGrid_.get_num_R_grid(); ++iR)
	{
		size_t dataGridIndex = spaceGrid_.R_conseq_local_to_data_conseq( iR );
		size_t thisFFTDataBegin = dataGridIndex*blockSize_*dimTimeFT_;
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,true);

		fftw_execute( fftw3PlanTimeBkwd_ );

		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,false);
	}
}

template<typename T>
void FFTBase<T>::perform_freq_to_time_fft_in_R()
{
	for ( size_t iR = 0 ; iR < spaceGrid_.get_num_R_grid(); ++iR)
	{
		size_t dataGridIndex = spaceGrid_.R_conseq_local_to_data_conseq( iR );
		size_t thisFFTDataBegin = dataGridIndex*blockSize_*dimTimeFT_;
		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,true);

		fftw_execute( fftw3PlanTimeFwd_ );

		this->copy_or_fill_buffer_time_FFT(thisFFTDataBegin,false);
	}
}

template<typename T>
void FFTBase<T>::copy_or_fill_buffer_first_FFT(size_t ikx, size_t indexTime, bool fillBuffer )
{
	auto FFT_layout = [&] ( size_t igpts, size_t ptsFFT, size_t ib )
	{
		T * bufferAsCmplxArray = reinterpret_cast<T*>(FFTBuffer_);
		T * thisLocation = bufferAsCmplxArray + ib*ptsFFT + igpts;
		return thisLocation;
	};

	size_t numFFTPts = spaceGrid_.get_num_k_grid()/spaceGrid_.get_k_grid().front();
	auto it = data_.begin();
	for ( size_t igpts = 0 ; igpts < numFFTPts ; ++igpts )
	{
		//We are in the k(R) grid layout and thus k(R)x is the slowest running dimension
		//and the current igpts index is related to the full k(R) index
		//ig by
		size_t ig = ikx*numFFTPts + igpts;
		size_t dataGridIndex = spaceGrid_.k_conseq_local_to_data_conseq(ig);
		size_t dataIndexBlock = (dataGridIndex*dimTimeFT_+indexTime)*blockSize_;

		for ( size_t ib = 0 ; ib < blockSize_ ; ++ib )
		{
			size_t index = dataIndexBlock+ib;
			if ( fillBuffer )
			{
				( *FFT_layout(igpts,numFFTPts,ib) ) = ( *(it+index) );
				//In debug mode, we check for NaNs. Before ...
				assert( (*(it+index)) == (*(it+index))  );
			}
			else
			{
				( *(it+index) ) = ( *FFT_layout(igpts,numFFTPts,ib) );
				//and after the FFT
				assert( (*(it+index)) == (*(it+index))  );
			}
		}
	}
}

template<typename T>
void FFTBase<T>::copy_or_fill_buffer_second_FFT(size_t iyz, size_t indexTime, bool fillBuffer )
{
	auto FFT_layout = [&] ( size_t igpts, size_t ptsFFT, size_t ib )
	{
		T * bufferAsCmplxArray = reinterpret_cast<T*>(FFTBuffer_);
		T * thisLocation = bufferAsCmplxArray+ib*ptsFFT + igpts;
		return thisLocation;
	};

	size_t numFFTPts = spaceGrid_.get_R_grid().front();

	auto it = data_.begin();
	for ( size_t igpts = 0 ; igpts < numFFTPts ; ++igpts )
	{
		//We are in the R grid layout and thus Rx is the fastest running dimension.
		//the current igpts index is related to the full R index
		//ig by
		size_t ig = iyz*numFFTPts + igpts;
		size_t dataRGridIndex = spaceGrid_.R_conseq_local_to_data_conseq(ig);
		//									^ see above
		size_t dataIndexBlock = (dataRGridIndex*dimTimeFT_+indexTime)*blockSize_;

		for ( size_t ib = 0 ; ib < blockSize_ ; ++ib )
		{
			size_t index = dataIndexBlock+ib;
			if ( fillBuffer )
			{
				( *FFT_layout(igpts,numFFTPts,ib) ) = ( *(it+index) );
				//In debug mode, we check for NaNs, Before
				assert( (*(it+index)) == (*(it+index))  );
			}
			else
			{
				( *(it+index) ) = ( *FFT_layout(igpts,numFFTPts,ib) );
				//... and after the FFT
				assert( (*(it+index)) == (*(it+index))  );
			}
		}
	}
}

template<typename T>
void FFTBase<T>::perform_k_to_R_fft()
{
	//First are the non-parallelized grid parts
	for ( size_t it = 0 ; it < this->get_num_time(); ++it)
	{
		//a one (or two) D transform for each kx point
		for (size_t ikx = 0; ikx < spaceGrid_.get_k_grid().front(); ++ikx)
		{
			this->copy_or_fill_buffer_first_FFT(ikx,it,true);
			fftw_execute( fftw3PlanGridParaBkwd_ );
			this->copy_or_fill_buffer_first_FFT(ikx,it,false);
		}
	}

	//Then we transpose on the grid
	spaceGrid_.grid_data_transposition(true,data_,dimTimeFT_*blockSize_);

	//and transform the parallelized grid parts
	for ( size_t it = 0 ; it < this->get_num_time(); ++it)
	{
		size_t dimRYZ = spaceGrid_.get_num_R_grid()/spaceGrid_.get_R_grid().front();
		//a one D transform along kx (still k!) for each Ry(and Rz) point
		for (size_t ikyz = 0; ikyz < dimRYZ; ++ikyz)
		{
			this->copy_or_fill_buffer_second_FFT(ikyz,it,true);
			fftw_execute( fftw3PlanGridSingleBkwd_ );
			this->copy_or_fill_buffer_second_FFT(ikyz,it,false);
		}
	}

	//Normalization
	for ( auto &d : data_ )
		d *= 1.0/spaceGrid_.get_num_grid();

	//In debug mode, we check for NaNs
	for ( auto d : data_ )
		assert( d == d );
}

template<typename T>
void FFTBase<T>::perform_R_to_k_fft()
{
	//First transform the non-parallelized grid parts. In this case along Rx.
	//We fill the buffer for each y (and z) with data for Rx and FFT
	for ( size_t it = 0 ; it < this->get_num_time(); ++it)
	{
		//a one D transform along Rx for each Ry(and Rz) point
		size_t dimYZ = spaceGrid_.get_num_R_grid()/spaceGrid_.get_R_grid().front();
		for (size_t iRx = 0; iRx < dimYZ; ++iRx)
		{
			this->copy_or_fill_buffer_second_FFT(iRx,it,true);
			fftw_execute( fftw3PlanGridSingleFwd_ );
			this->copy_or_fill_buffer_second_FFT(iRx,it,false);
		}
	}

	//Then we transpose on the grid (from R to k)
	spaceGrid_.grid_data_transposition(false,data_,dimTimeFT_*blockSize_);

	//and transform the parallelized grid parts
	for ( size_t it = 0 ; it < this->get_num_time(); ++it)
	{
		//a one (or two) D transform for each kx point
		for (size_t ikx = 0; ikx < spaceGrid_.get_k_grid().front(); ++ikx)
		{
			this->copy_or_fill_buffer_first_FFT(ikx,it,true);
			fftw_execute( fftw3PlanGridParaBkwd_ );
			this->copy_or_fill_buffer_first_FFT(ikx,it,false);
		}
	}
}

template<typename T>
T const * FFTBase<T>::read_phs_grid_ptr_block(size_t iGrid, size_t it ) const
{
	size_t igdata;
	if ( isKSpace_ )
		igdata = spaceGrid_.k_conseq_local_to_data_conseq( iGrid );
	if ( ! isKSpace_ )
		igdata = spaceGrid_.R_conseq_local_to_data_conseq( iGrid );
	size_t dIndex = (igdata*dimTimeFT_+it)*blockSize_;
	assert( dIndex < data_.size() );
	return &(data_[ dIndex ]);
}

template<typename T>
T * FFTBase<T>::write_phs_grid_ptr_block(size_t iGrid, size_t it )
{
	size_t igdata;
	if ( isKSpace_ )
		igdata = spaceGrid_.k_conseq_local_to_data_conseq( iGrid );
	if ( ! isKSpace_ )
		igdata = spaceGrid_.R_conseq_local_to_data_conseq( iGrid );
	size_t dIndex = (igdata*dimTimeFT_+it)*blockSize_;
	assert( dIndex < data_.size() );
	return &(data_[ dIndex ]);
}

template<typename T>
T const * FFTBase<T>::read_data_ptr_block(size_t iGrid, size_t it ) const
{
	return &(data_[ (iGrid*dimTimeFT_+it)*blockSize_ ]);
}

template<typename T>
T * FFTBase<T>::write_data_ptr_block(size_t iGrid, size_t it )
{
	return &(data_[ (iGrid*dimTimeFT_+it)*blockSize_ ]);
}

template<typename T>
bool FFTBase<T>::is_init() const
{
	return isInit_;
}

template<typename T>
bool FFTBase<T>::is_in_k_space() const
{
	return isKSpace_;
}

template<typename T>
bool FFTBase<T>::is_in_time_space() const
{
	return isTimeSpace_;
}

template<typename T>
void FFTBase<T>::set_time_space(bool newState)
{
	isTimeSpace_ = newState;
}

template<typename T>
void FFTBase<T>::perform_time_fft()
{
	if ( isTimeSpace_ )
	{
		if ( isKSpace_ )
			this->perform_time_to_freq_fft_in_k();

		if ( ! isKSpace_ )
			this->perform_time_to_freq_fft_in_R();
	}
	if ( ! isTimeSpace_ )
	{
		if ( isKSpace_ )
			this->perform_freq_to_time_fft_in_k();

		if ( ! isKSpace_ )
			this->perform_freq_to_time_fft_in_R();
	}
	isTimeSpace_ = ! isTimeSpace_;
}

template<typename T>
void FFTBase<T>::perform_space_fft()
{
	if ( isKSpace_ )
		this->perform_k_to_R_fft();

	if ( ! isKSpace_ )
		this->perform_R_to_k_fft();

	isKSpace_ = ! isKSpace_;
}

} /* namespace gw_flex */
} /* namespace scallop */
