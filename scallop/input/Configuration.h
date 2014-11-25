/*	This file Configuration.h is part of scallop.
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

#ifndef SCALLOP_INPUT_CONFIGURATION_H_
#define SCALLOP_INPUT_CONFIGURATION_H_

#include "scallop/input/InputBase.h"

namespace scallop {
namespace input {

class Configuration : public InputBase<Configuration> {

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			method,
			"Decide which method is used to compute SC\n"
			"Possible choices are Elishberg 'Eliash'\n"
			"\tand (Spin)SCDFT 'SCDFT'",
			"SCDFT",
			"SCDFT",
			std::string);
};

} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_CONFIGURATION_H_ */
