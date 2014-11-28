/*	This file EliashbergImaginaryAxis.h is part of scallop.
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

#ifndef SCALLOP_ELIASHBERG_ELIASHBERGIMAGINARYAXIS_H_
#define SCALLOP_ELIASHBERG_ELIASHBERGIMAGINARYAXIS_H_

#include "scallop/eliashberg/EliashbergGapFunction.h"
#include "scallop/eliashberg/NambuOffDiagonalEnergyIntegralN.h"
#include "scallop/eliashberg/FrequencyCorrectionZ.h"
#include "scallop/eliashberg/GapSquareRoot.h"

namespace scallop {
namespace eliashberg {

/**
 * \brief Solve the Eliashberg equations at given temperature and splitting
 */
template<typename T>
class EliashbergImaginaryAxis {

	void solve(
			T _temperature,
			EliashbergGapFunction<T> const& deltaEInitalGuess,
			FrequencyCorrectionZ<T> const& ZInitalGuess);

private:

	T _temperature;

	EliashbergGapFunction<T> _deltaE;
	FrequencyCorrectionZ<T> _Z;
	GapSquareRoot<T> _gapSquareRootSpinUp;
	GapSquareRoot<T> _gapSquareRootSpinDown;


};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/EliashbergImaginaryAxis.hpp"
#endif /* SCALLOP_ELIASHBERG_ELIASHBERGIMAGINARYAXIS_H_ */
