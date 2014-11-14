/*	This file Setup.h is part of scallop.
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
 *  Created on: Nov 14, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_INPUT_SETUP_H_
#define SCALLOP_INPUT_SETUP_H_

namespace scallop {
namespace input {

/**
 * Perform the initial startup of the code.
 */
class Setup {
public:

	/**
	 * Create the initial startup.
	 *
	 * @param argc The number of arguments taken from the input.
	 * @param argv The array of size argc with c-stype char arrays.
	 */
	Setup(int argc, char *argv[]);
private:
};

} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_SETUP_H_ */
