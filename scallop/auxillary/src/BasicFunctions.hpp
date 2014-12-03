/*	This file BasicFunctions.hpp is part of scallop.
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

#include "scallop/auxillary/BasicFunctions.h"
#include <cmath>

namespace scallop {
namespace auxillary {

template<typename T>
T BasicFunctions::inverse_temperature(T temperature){
	return 1/(Constants<T>::kBoltzmannHartree*temperature);
}

template<typename T>
T BasicFunctions::matzubara_frequency_of_index(int index, T beta) {
	return M_PI*(2*index+1)/beta;
}

} /* namespace auxillary */
} /* namespace scallop */