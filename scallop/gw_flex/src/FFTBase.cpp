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
	//We assume that the buffer were properly filled with the time and block data for each grid point
	int * n = new int [1];
	n[0] = dimTimeFT_;
	int * inembed =NULL;
	int * onembed = NULL;
	int howmany = static_cast<int>( blockSize_ );
	int idist = static_cast<int>(dimTimeFT_);
	int odist = static_cast<int>(dimTimeFT_);
	int istride = 1;
	int ostride = 1;

	fftw3PlanTimeFwd_ = fftw_plan_many_dft(
			1,n,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			-1,FFTW_PATIENT);

	fftw3PlanTimeBkwd_ = fftw_plan_many_dft(
			1,n,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			1,FFTW_PATIENT);

	delete [] n;
	delete [] inembed;
	delete [] onembed;
}


template<>
void
FFTBase< std::complex<double> >::plan_space_fft()
{
	//We do the first k to R space transform; we do all FFT for the blockSize_ at once
	//We assume that the buffer is properly filled with the non-parallel space dims and block data for each grid point
	//Note that k is parallelized in the k direction. Thus our first space FFT is the y (and possibly z)
	//The assumed layout in the buffer is blocksize_ [space dim-1] dimensional block of data values for a given block and time index
	int dim = spaceGrid_.get_dim()-1;
	int * n = new int [ dim ];
	for ( int id = 0 ; id < dim; ++id)
		n[id] = spaceGrid_.get_k_grid()[id+1];
	int dblock = 1;
	for ( int id = 0 ; id < dim; ++id)
		dblock*=n[id];
	int * inembed =NULL;
	int * onembed = NULL;
	int howmany = static_cast<int>( blockSize_ );
	int idist = dblock;
	int odist = dblock;
	int istride = 1;
	int ostride = 1;

	fftw3PlanGridParaBkwd_ = fftw_plan_many_dft(
			dim,n,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			1,FFTW_PATIENT);

	//Check the second k to R space transform; we do all FFT for the blockSize_ at once
	//Our second FFT is the x direction
	//The assumed layout in the buffer is filled with
	// blocksize_ [first total space dim] dimensional block of data values for a given block and time index
	int n1d = spaceGrid_.get_grid().front();
	idist = odist = n1d;
	fftw3PlanGridSingleBkwd_ = fftw_plan_many_dft(
			1,&n1d,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			1,FFTW_PATIENT);

	//We plan the first R to k space transform; we do all FFT for the blockSize_ at once
	//This plan is for the Rx dimension
	fftw3PlanGridSingleFwd_ = fftw_plan_many_dft(
			1,&n1d,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			-1,FFTW_PATIENT);

	//We plan the second R to k space transform; we do all FFT for the blockSize_ at once
	//Note that R is parallelized in the last direction.
	idist = odist = dblock;
	fftw3PlanGridParaFwd_ = fftw_plan_many_dft(
			dim,n,howmany,
			FFTBuffer_, inembed,istride, idist,
			FFTBuffer_, onembed,ostride, odist,
			-1,FFTW_PATIENT);


	delete [] n;
	delete [] inembed;
	delete [] onembed;
}

} /* namespace gw_flex */
} /* namespace scallop */
