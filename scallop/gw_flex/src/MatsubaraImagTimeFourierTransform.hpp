/*	This file MatsubaraImagTimeFourierTransform.hpp is part of scallop.
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
 *  Created on: Oct 28, 2016
 *      Author: Andreas Linscheid
 */

#include "scallop/gw_flex/MatsubaraImagTimeFourierTransform.h"
#include "scallop/error_handling/error_handling.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
MatsubaraImagTimeFourierTransform<T>::MatsubaraImagTimeFourierTransform( bool Fermi )
	: isFermi_(Fermi)
{

};

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::initialize(
		size_t dimImTime,
		std::vector<size_t> gridDims,
		size_t dataBlockSize,
		bool initialInFreqDomain,
		bool initialInReciprocalDomain,
		std::vector<T> data)
{
	size_t nK = 1;
	for (auto kd : gridDims )
		nK *= kd;
	FFTBase<T>::plan_ffts(
			std::move(data),
			gridDims,
			dimImTime,
			dataBlockSize);
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::transform_Mfreq_to_itime( bT invTemp )
{
	if ( statusTimeDomain_ )
		scallop::error_handling::Error(
				"Cannot apply frequency to time transform to something already in time space!");

	if ( isFermi_ )
	{
		this->fourier_transform_fermions_freq_to_time( invTemp );
	}
	else
	{
		this->fourier_transform_bosons_freq_to_time( invTemp );
	}
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::transform_itime_to_Mfreq( bT invTemp )
{
	if ( ! statusTimeDomain_ )
		scallop::error_handling::Error(
				"Cannot apply time to frequency transform to something already in frequency space!");

	if ( isFermi_ )
	{
		this->fourier_transform_fermions_time_to_freq(invTemp);
	}
	else
	{
		this->fourier_transform_bosons_time_to_freq(invTemp);
	}
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_fermions_freq_to_time( bT invTemp )
{
	//frequency to time transform is a simple DFT, normalized by the inverse temperature
	this->perform_freq_to_time_fft();

	//Multiply the phase factor exp( - \sqrt(-1) \pi i/N_mats ) from the (2n+1)
	for (size_t ik = 0 ; ik < this->get_num_grid() ; ++ik )
		for (size_t iw = 0 ; iw < this->get_num_time() ; ++iw )
		{
			auto it = this->data_begin_modify() + (ik*this->get_num_time()+iw)*this->get_data_block_size();
			auto itEnd = it+this->get_data_block_size();
			for ( ; it != itEnd ; ++it )
				(*it) *= std::exp( T(0, -(M_PI * iw) / this->get_num_time() ) );
		}

	//Scale by beta
	for ( auto it = this->data_begin_modify(); it != this->data_end_modify(); ++it  )
		(*it) *= 1.0/invTemp;

	statusTimeDomain_ = true;
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_bosons_freq_to_time( bT invTemp )
{

}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_fermions_time_to_freq( bT invTemp )
{
	const size_t nM = this->get_num_time();

	//Multiply the phase factor exp( \sqrt(-1) \pi i/N_mats ) from the (2n+1)
	for (size_t ik = 0 ; ik < this->get_num_grid() ; ++ik )
		for (size_t iw = 0 ; iw < nM ; ++iw )
		{
			auto it = this->data_begin_modify() + (ik*nM+iw)*this->get_data_block_size();
			auto itEnd = it+this->get_data_block_size();
			for ( ; it != itEnd ; ++it )
				(*it) *=  std::exp( T(0, (M_PI * iw) / nM) );
		}

	this->perform_time_to_freq_fft();

	//Multiply prefactor that comes from the integration kernel
	for (size_t ik = 0 ; ik < this->get_num_grid() ; ++ik )
		for (size_t iw = 0 ; iw < this->get_num_time() ; ++iw )
		{
			auto it = this->data_begin_modify() + (ik*nM+iw)*this->get_data_block_size();
			auto itEnd = it+this->get_data_block_size();

			//remember that frequencies are stored with negative frequencies in the second half of the array
			int frequencyIndex =
					(iw < nM/2 ?
							static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

			T phase = T(0,M_PI*(2*frequencyIndex+1));
			for ( ; it != itEnd ; ++it )
				(*it) *=  -invTemp*(1.0-std::exp(phase/bT(nM)))/phase;
		}
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_bosons_time_to_freq( bT invTemp )
{

}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::set_time_domain()
{
	statusTimeDomain_ = true;
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::set_frequency_domain()
{
	statusTimeDomain_ = false;
}

} /* namespace gw_flex */
} /* namespace scallop */
