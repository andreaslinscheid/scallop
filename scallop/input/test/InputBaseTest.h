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

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(sizetTest,
			"Tests the reading of a type of size_t.\n",
			zero,0,size_t);

//	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(boolTest,
//			"Tests the reading of a type of type bool.\n"
//			"The default is false",
//			false,bool);
//
//	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(doubleTest,
//			"Tests the reading of a type of type double.\n"
//			"The default is 4.0 .",
//			4.0,double);
//
//
//	std::vector<size_t> defaultSizeTVect{10,60,30};
//
//	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(vectorSizeTTest,
//			"Tests the reading of a vector of type size_t.\n"
//			"The default is 10,60,30",
//			,
//			std::vector<size_t>);
private:
//	//test several types of input
//	size_t stypeOption;
//	bool _yesNoOption;
//	double _floatDOption;
//	float _floatOption;
//	std::vector<size_t> _vectorOfNumbersOption;
//	std::string _textOption;
};

} /* namespace test */
} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_TEST_INPUTBASETEST_H_ */
