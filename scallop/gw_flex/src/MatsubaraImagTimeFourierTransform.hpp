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
		bool initialInTimeDomain,
		bool initialInReciprocalDomain,
		typename auxillary::TemplateTypedefs<T>::scallop_vector data)
{
	FFTBase<T>::initialize(
			std::move(data),
			initialInReciprocalDomain,
			initialInTimeDomain,
			gridDims,
			dimImTime,
			dataBlockSize);
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::transform_itime_Mfreq( bT invTemp )
{
	if ( isFermi_ )
	{
		this->fourier_transform_fermions_time_freq( invTemp );
	}
	else
	{
		this->fourier_transform_bosons_time_freq( invTemp );
	}
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_fermions_time_freq( bT invTemp )
{
	const size_t nM = this->get_num_time();
	const size_t nK = this->get_spaceGrid_proc().get_num_k_grid();
//	const size_t nR = this->get_spaceGrid_proc().get_num_R_grid();

	if ( this->is_in_time_space() )
	{
		if ( this->is_in_k_space() )
		{
			//Multiply the phase factor exp( \sqrt(-1) \pi i/N_mats ) from the (2n+1)
			for (size_t ik = 0 ; ik < nK ; ++ik )
			{
				size_t id = this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ik);
				for (size_t iw = 0 ; iw < nM ; ++iw )
				{
					auto ptr_block = this->write_data_ptr_block(id,iw);
					for ( size_t ib = 0; ib < this->get_data_block_size(); ++ib)
						ptr_block[ib] *=  std::exp( T(0, (M_PI * iw) / nM) );
				}
			}
		}
	}

	//frequency to time transform is a simple DFT, normalized by the inverse temperature
	this->perform_time_fft();

	if ( ! this->is_in_time_space() )
	{
		if ( this->is_in_k_space() )
		{
			//Multiply prefactor that comes from the integration kernel
			for (size_t ik = 0 ; ik < nK ; ++ik )
			{
				size_t id = this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ik);
				for (size_t iw = 0 ; iw < nM ; ++iw )
				{
					auto ptr_block = this->write_data_ptr_block(id,iw);

					//remember that frequencies are stored with negative frequencies in the second half of the array
					int frequencyIndex =
							(iw < nM/2 ?
									static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

					T phase = T(0,M_PI*(2*frequencyIndex+1));
					for ( size_t ib = 0; ib < this->get_data_block_size(); ++ib)
						ptr_block[ib] *=  -invTemp*(1.0-std::exp(phase/bT(nM)))/phase;
				}
			}
		}
	}
	else //... we are now in time space
	{
		//Multiply the phase factor exp( - \sqrt(-1) \pi i/N_mats ) from the (2n+1)
		// and scale by beta
		if ( this->is_in_k_space() )
		{
			for (size_t ik = 0 ; ik < nK ; ++ik )
			{
				size_t idk = this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ik);
				for (size_t iw = 0 ; iw < nM ; ++iw )
				{
					auto ptr_block = this->write_data_ptr_block(idk,iw);
					for ( size_t ib = 0; ib < this->get_data_block_size(); ++ib)
						ptr_block[ib] *= std::exp( T(0, -(M_PI * iw) / nM ) )/invTemp;
				}
			}
		}
	}

}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::fourier_transform_bosons_time_freq( bT invTemp )
{

}

} /* namespace gw_flex */
} /* namespace scallop */
