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

	this->initialize_corrections();
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::initialize_corrections()
{
	typedef std::complex<double> CD;

	size_t nM = this->get_num_time();

	aN1_ = aN0_ = bN1_ = bN0_ = WN_ =
			typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector(nM);

	for ( size_t iw = 0; iw < nM; ++iw )
	{
		int frequencyIndex =
				(iw < nM/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );
		double pn = M_PI*( isFermi_ ? 2*frequencyIndex + 1 : 2 * frequencyIndex  ) / double(nM);

		//See if we use the Taylor series
		if ( std::abs(pn) < 0.2 )
		{
			WN_[iw] = 1 - std::pow(pn,4)/60. + (13*std::pow(pn,6))/7560.;

			aN0_[iw] = 0.041666666666666664 - CD(0,0.17447916666666666)*pn + std::pow(pn,2)/8. + CD(0,0.0886501736111111)*std::pow(pn,3)
						- (11089*std::pow(pn,4))/322560. - CD(0,0.014633614676339286)*std::pow(pn,5) + (214993*std::pow(pn,6))/4.644864e7 +
					   CD(0,0.0014226102744881227)*std::pow(pn,7);

			aN1_[iw] = 0.041666666666666664 + CD(0,0.453125)*pn - (193*std::pow(pn,2))/384. - CD(0,0.3869357638888889)*std::pow(pn,3)
						+ (25987*std::pow(pn,4))/107520. + CD(0,0.12823331318204365)*std::pow(pn,5) - (787147*std::pow(pn,6))/1.327104e7
						- CD(0,0.0243327226588335)*std::pow(pn,7);


			bN1_[iw] = 0.041666666666666664 - CD(0,0.453125)*pn - (193*std::pow(pn,2))/384. + CD(0,0.3869357638888889)*std::pow(pn,3)
						+ (25987*std::pow(pn,4))/107520. - CD(0,0.12823331318204365)*std::pow(pn,5) - (787147*std::pow(pn,6))/1.327104e7
						+ CD(0,0.0243327226588335)*std::pow(pn,7);

			bN0_[iw] =0.041666666666666664 + CD(0,0.17447916666666666)*pn + std::pow(pn,2)/8. - CD(0,0.0886501736111111)*std::pow(pn,3)
						- (11089*std::pow(pn,4))/322560. + CD(0,0.014633614676339286)*std::pow(pn,5) + (214993*std::pow(pn,6))/4.644864e7 -
					   CD(0,0.0014226102744881227)*std::pow(pn,7);
		}
		else // no we don't use the series
		{
			WN_[iw]  = (4*std::pow(std::sin(pn/2.),3)*(2*std::cos(pn/2.) + pn*std::sin(pn/2.)))/std::pow(pn,3);

			aN0_[iw] = (CD(0,4.0) - 2.0*pn - 12.0*std::exp(CD(0,2.0*pn))*pn + 8.0*std::exp(CD(0,pn))*(CD(0,-1.0) + pn)
						+ 8.0*std::exp(CD(0,3.0*pn))*(CD(0,1.0) + pn)
						+ 2.0*std::exp(CD(0,4.0*pn))*(CD(0,2.0) + pn) + std::exp(CD(0,(3.0*pn)/2.0))*(CD(0,-8.0) + (16.0 + CD(0,15.0)*pn)*pn))/
					   (8.0*std::exp(CD(0,(3.0*pn)/2.0))*std::pow(pn,3));

			aN1_[iw] = -(CD(0,-2.0) + pn + 6.0*std::exp(CD(0,2.0*pn))*pn - 4.0*std::exp(CD(0,pn))*(CD(0,-1.0) + pn)
						+ 4.0*std::exp(CD(0,3.0*pn))*(CD(0,1.0) + pn) + std::exp(CD(0,4.0*pn))*(CD(0,2.0) + pn)
						+std::exp(CD(0,pn/2.0))*(CD(0,-8.0) + (12.0 + CD(0,5.0)*pn)*pn))/(4.0*std::exp(CD(0,pn/2.0))*std::pow(pn,3));

			bN1_[iw] = -(CD(0,-2.0) + pn + 6.0*std::exp(CD(0,2.0*pn))*pn + 4.0*std::exp(CD(0,pn))*(CD(0,-1.0) + pn)
						- 4.0*std::exp(CD(0,3*pn))*(CD(0,1.0) + pn) + std::exp(CD(0,4.0*pn))*(CD(0,2.0) + pn)
						+ std::exp(CD(0,7.0*pn/2.0))*(CD(0,8.0) + (12.0 - CD(0,5.0)*pn)*pn))
						/(4.0*std::exp(CD(0,7.0*pn/2.0))*std::pow(pn,3));

			bN0_[iw] = (-12.0*std::exp(CD(0,2.0*pn))*pn + 8.0*std::exp(CD(0,pn))*(CD(0,-1.0) + pn) + 8.0*std::exp(CD(0,3.0*pn))*(CD(0,1.0) + pn)
						+ 2.0*(CD(0,-2.0) + pn) - 2.0*std::exp(CD(0,4.0*pn))*(CD(0,2.0) + pn) + std::exp(CD(0,5.0*pn/2.0))*(CD(0,8.0)
						+ (16.0 - CD(0,15.0)*pn)*pn))/(8.0*std::exp(CD(0,5.0*pn/2.0))*std::pow(pn,3));

		}
		bN1_[iw] *= std::exp( CD(0,pn*nM) );
		bN0_[iw] *= std::exp( CD(0,pn*nM) );
	}

	size_t nGR = this->get_spaceGrid_proc().get_num_R_grid();
	size_t nGK = this->get_spaceGrid_proc().get_num_k_grid();
	size_t nGm = nGK > nGR ? nGK : nGR;
	buffFNMm2_ = buffFNMm1_ = buffF1_ = buffF0_ = std::vector<T>(nGm*this->get_data_block_size());
}

template<typename T>
void MatsubaraImagTimeFourierTransform<T>::transform_itime_Mfreq( bT invTemp )
{
	const size_t nM = this->get_num_time();
	const size_t nG = this->is_in_k_space() ?
			this->get_spaceGrid_proc().get_num_k_grid() :
			this->get_spaceGrid_proc().get_num_R_grid() ;
	bool timeToFreqTransform = this->is_in_time_space();
	const size_t nB = this->get_data_block_size();

	if ( timeToFreqTransform )
	{
		//Multiply the phase factor exp( \sqrt(-1) \pi i/N_mats ) from the (2n+1)
		for (size_t ig = 0 ; ig < nG ; ++ig )
		{
			size_t id = this->is_in_k_space() ?
					this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ig) :
					this->get_spaceGrid_proc().R_conseq_local_to_data_conseq(ig) ;

			if ( isFermi_ )
			{
				for (size_t it = 0 ; it < nM ; ++it )
				{
					auto ptr_block = this->write_data_ptr_block(id,it);
					for ( size_t ib = 0; ib < nB; ++ib)
						ptr_block[ib] *=  std::exp( T(0, (M_PI * it) / nM) );
				}
			}

			for ( size_t ib = 0; ib < nB; ++ib)
			{
				//Note: we have to ensure nM >= 2
				assert( nM >= 2 );
				buffF0_[ig*nB+ib] = this->read_data_ptr_block(id,0)[ib];
				buffF1_[ig*nB+ib] = this->read_data_ptr_block(id,1)[ib];;
				buffFNMm1_[ig*nB+ib] = this->read_data_ptr_block(id,nM-1)[ib];;
				buffFNMm2_[ig*nB+ib] = this->read_data_ptr_block(id,nM-2)[ib];;
			}
		}
	}
	else //We transform from frequency to time
	{
		//Multiply the phase factor exp( \sqrt(-1) \pi n/N_mats )
		for (size_t ig = 0 ; ig < nG ; ++ig )
		{
			size_t id = this->is_in_k_space() ?
					this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ig) :
					this->get_spaceGrid_proc().R_conseq_local_to_data_conseq(ig) ;

			for (size_t iw = 0 ; iw < nM ; ++iw )
			{
				int frequencyIndex =
						(iw < nM/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

				auto ptr_block = this->write_data_ptr_block(id,iw);
				for ( size_t ib = 0; ib < nB; ++ib)
					ptr_block[ib] *=  std::exp( T(0, -(M_PI * frequencyIndex) / nM) );
			}
		}
	}

	//frequency to time transform is a simple DFT, normalized by the inverse temperature
	this->perform_time_fft();

	//Multiply prefactor that comes from the integration kernel and the end-point corrections
	if ( timeToFreqTransform )
	{
		for (size_t ig = 0 ; ig < nG ; ++ig )
		{
			size_t id = this->is_in_k_space() ?
					this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ig) :
					this->get_spaceGrid_proc().R_conseq_local_to_data_conseq(ig) ;

			for (size_t iw = 0 ; iw < nM ; ++iw )
			{
				int frequencyIndex =
						(iw < nM/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

				auto ptr_block = this->write_data_ptr_block(id,iw);
				for ( size_t ib = 0; ib < nB; ++ib)
					ptr_block[ib] =  (invTemp/nM) * (ptr_block[ib]*WN_[iw]*
							std::exp( T(0, (M_PI * (isFermi_ ? frequencyIndex+0.5 : frequencyIndex) ) / nM) ) );
//					    - aN0_[iw]*buffF0_[ig*nB+ib] - aN1_[iw]*buffF1_[ig*nB+ib]
//						- bN1_[iw]*buffFNMm2_[ig*nB+ib] - bN0_[iw]*buffFNMm1_[ig*nB+ib] );
			}
		}
	}
	else //timeToFreqTransform is false and we transform frequency to time
	{
		//Multiply the phase factor exp( - \sqrt(-1) \pi (i+0.5)/N_mats ) from the (2n+1)
		// and scale by beta
		for (size_t ig = 0 ; ig < nG ; ++ig )
		{
			size_t id = this->is_in_k_space()  ?
					this->get_spaceGrid_proc().k_conseq_local_to_data_conseq(ig) :
					this->get_spaceGrid_proc().R_conseq_local_to_data_conseq(ig) ;
			for (size_t it = 0 ; it < nM ; ++it )
			{
				auto ptr_block = this->write_data_ptr_block(id,it);
				for ( size_t ib = 0; ib < nB; ++ib)
					ptr_block[ib] *= isFermi_ ?
							std::exp( T(0, -(M_PI * (it+0.5)) / nM ) )/invTemp :
							1.0/invTemp;
			}
		}
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
