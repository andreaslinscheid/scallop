/*	This file MatsubaraImagTimeFourierTransform.h is part of scallop.
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
 *      Author: alinsch
 */

#ifndef SCALLOP_GW_FLEX_MATSUBARAIMAGTIMEFOURIERTRANSFORM_H_
#define SCALLOP_GW_FLEX_MATSUBARAIMAGTIMEFOURIERTRANSFORM_H_

#include "scallop/gw_flex/FFTBase.h"
#include "scallop/auxillary/TypeMapComplex.h"
#include <vector>
#include <cstdlib>

namespace scallop {
namespace gw_flex {

/**
 *	Module to transform data from the imaginary time domain to the
 *	Matsubara frequency domain and back.
 */
template<typename T>
class MatsubaraImagTimeFourierTransform : public FFTBase<T>
{
	typedef typename scallop::auxillary::TypeMapComplex<T>::type bT;

public:

	MatsubaraImagTimeFourierTransform( bool Fermi );

	void initialize(
			size_t dimImTime,
			std::vector<size_t> gridDims,
			size_t dataBlockSize,
			bool initialInFreqDomain,
			bool initialInReciprocalDomain,
			std::vector<T> data);

	/**
	 * Transform \p data from imaginary time to Matsubara frequency.
	 *
	 * @param data	Data to be operated on. We assume a layout [\p howmany][it][\p blockSize], such that
	 * 				two consecutive data elements belonging to the same entity in time are \p blockSize apart.
	 * @param howmany
	 * @param blockSize
	 * @param Fermi
	 */
	void transform_itime_to_Mfreq( bT invTemp );

	void transform_Mfreq_to_itime( bT invTemp );

	/**
	 * 	When called overwrites the internal status such that
	 * 	the data is interpreted to be in the time domain
	 */
	void set_time_domain();

	/**
	 * 	When called overwrites the internal status such that
	 * 	the data is interpreted to be in the frequency domain
	 */
	void set_frequency_domain();

private:

	//If set true, it will call the time Fourier transform for Fermions, otherwise the one for Bosons
	bool isFermi_;

	bool isInit_ = false;

	//track in which domain this object currently resides
	bool statusTimeDomain_ = false;

	void fourier_transform_fermions_time_to_freq( bT invTemp );

	void fourier_transform_fermions_freq_to_time( bT invTemp );

	void fourier_transform_bosons_time_to_freq( bT invTemp );

	void fourier_transform_bosons_freq_to_time( bT invTemp );


};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/MatsubaraImagTimeFourierTransform.hpp"
#endif /* SCALLOP_GW_FLEX_MATSUBARAIMAGTIMEFOURIERTRANSFORM_H_ */
