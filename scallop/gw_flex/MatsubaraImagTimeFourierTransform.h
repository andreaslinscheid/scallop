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
#include "scallop/auxillary/TemplateTypedefs.h"
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
			bool initialInTimeDomain,
			bool initialInReciprocalDomain,
			typename auxillary::TemplateTypedefs<T>::scallop_vector data);

	/**
	 * Transform \p data from imaginary time to Matsubara frequency.
	 *
	 * @param data	Data to be operated on. We assume a layout [\p howmany][it][\p blockSize], such that
	 * 				two consecutive data elements belonging to the same entity in time are \p blockSize apart.
	 * @param howmany
	 * @param blockSize
	 * @param Fermi
	 */
	void transform_itime_Mfreq( bT invTemp );

private:

	///If set true, it will call the time Fourier transform for Fermions, otherwise the one for Bosons
	bool isFermi_;

	///Buffer for the weight factors W(n) for the transform it->w_n
	///independent on T, this needs to be double at least
	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector WN_;

	///Endpoint corrections for the transform it->w_n
	///independent on T, this needs to be double at least
	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector aN0_;
	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector aN1_;
	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector bN0_;
	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector bN1_;

	std::vector<T> buffF0_;
	std::vector<T> buffF1_;
	std::vector<T> buffFNMm1_;
	std::vector<T> buffFNMm2_;

	void initialize_corrections();
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/MatsubaraImagTimeFourierTransform.hpp"
#endif /* SCALLOP_GW_FLEX_MATSUBARAIMAGTIMEFOURIERTRANSFORM_H_ */
