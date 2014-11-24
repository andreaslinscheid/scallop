/*	This file InputBaseTest.h is part of scallop.
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
 *  Created on: Nov 15, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_INPUT_TEST_INPUTBASETEST_H_
#define SCALLOP_INPUT_TEST_INPUTBASETEST_H_

#include "scallop/input/InputBase.h"
#include <vector>

namespace scallop {
namespace input {
namespace test {

class InputBaseTest : public InputBase<InputBaseTest> {

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			sizeTest,
			"Tests the reading of a type of size_t.\n",
			"zero",
			0,
			size_t);

	INPUTBASE_INPUT_OPTION_MACRO(
			doubleTest,
			"Tests the reading of a type of double.\n",
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			vectorSizeT,
			"Tests the reading of a vector of size_ts.\n",
			"List of 0 11 and 42",
			{0 COMMA_SUBSTITUTION 11 COMMA_SUBSTITUTION 42},
			std::vector<size_t>);

	INPUTBASE_INPUT_OPTION_MACRO(
			boolTest,
			"Tests the bool reading\n",
			bool);
};

} /* namespace test */
} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_TEST_INPUTBASETEST_H_ */
