/*	This file FrequencyCorrectionZ.h is part of scallop.
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
 *  Created on: Nov 27, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_FREQUENCYCORRECTIONZ_H_
#define SCALLOP_ELIASHBERG_FREQUENCYCORRECTIONZ_H_

namespace scallop {
namespace eliashberg {

template<typename T>
class FrequencyCorrectionZ : public MatzubaraSplittingVector<T> {

	void compute(
			MatzubaraEffectiveCouplingMatrix<T> const& EffectiveCouplingKernelUpSpin,
			MatzubaraEffectiveCouplingMatrix<T> const& EffectiveCouplingKernelDownSpin,
			NambuDiagonalEnergyIntegralM<T,auxillary::Constants<T>::upspin> const& MUpSpin,
			NambuDiagonalEnergyIntegralM<T,auxillary::Constants<T>::downspin> const& MDownSpin,
			FrequencyCorrectionZ<T> const& previousZ,
			double inverseTemperature);
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/FrequencyCorrectionZ.hpp"
#endif /* SCALLOP_ELIASHBERG_FREQUENCYCORRECTIONZ_H_ */
