/*	This file MatzubaraEffectiveCouplingMatrix.h is part of scallop.
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
 *  Created on: Nov 25, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_MATZUBARAEFFECTIVECOUPLINGMATRIX_H_
#define SCALLOP_ELIASHBERG_MATZUBARAEFFECTIVECOUPLINGMATRIX_H_

#include "scallop/eliashberg/MatzubaraSplittingMatrix.h"
#include "scallop/couplings/FrequencyDependentCoupling.h"

namespace scallop {
namespace eliashberg {

template<typename T>
class MatzubaraEffectiveCouplingMatrix : public MatzubaraSplittingVector<T> {
public:

	T & operator() (size_t b, size_t j, size_t bp, size_t jp, int nMinusNp);

	T operator() (size_t b, size_t j, size_t bp, size_t jp, int nMinusNp) const;

	void compute(
			size_t numberMatzubaraPts,
			T temperature,
			couplings::FrequencyDependentCoupling<T> const& couplingData);
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/MatzubaraEffectiveCouplingMatrix.hpp"
#endif /* SCALLOP_ELIASHBERG_MATZUBARAEFFECTIVECOUPLINGMATRIX_H_ */
