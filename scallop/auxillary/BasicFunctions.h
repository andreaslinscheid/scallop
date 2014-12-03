/*	This file BasicFunctions.h is part of scallop.
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

#ifndef SCALLOP_AUXILLARY_BASICFUNCTIONS_H_
#define SCALLOP_AUXILLARY_BASICFUNCTIONS_H_

namespace scallop {
namespace auxillary {

struct BasicFunctions {

	template<typename T>
	static T inverse_temperature(T temperature);

	template<typename T>
	static T matzubara_frequency_of_index(int index, T beta);
};

template<typename T>
struct Constants {

	static constexp int upspin = 1;

	static const int downspin = -1;

	static const T kBoltzmannMHartree = 0.003166811429;

	static const T THZToMHartree = 0.151983;

	static const T RYToMHartree = 500;

	static const T cmToTheMinus1ToMHartree = 0.0000045563352812*1000;

};

} /* namespace auxillary */
} /* namespace scallop */
#include "scallop/auxillary/src/BasicFunctions.hpp"
#endif /* SCALLOP_AUXILLARY_BASICFUNCTIONS_H_ */