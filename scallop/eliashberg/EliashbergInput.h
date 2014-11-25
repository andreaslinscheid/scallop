/*	This file EliashbergInput.h is part of scallop.
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

#ifndef SCALLOP_ELIASHBERG_ELIASHBERGINPUT_H_
#define SCALLOP_ELIASHBERG_ELIASHBERGINPUT_H_

#include "scallop/input/InputBase.h"

namespace scallop {
namespace eliashberg {

class EliashbergInput : public input::InputBase<EliashbergInput> {

	INPUTBASE_INPUT_OPTION_MACRO(temp,
			"The temperature of the system that we calculate",
			double);
};

} /* namespace eliashberg */
} /* namespace scallop */
#endif /* SCALLOP_ELIASHBERG_ELIASHBERGINPUT_H_ */
