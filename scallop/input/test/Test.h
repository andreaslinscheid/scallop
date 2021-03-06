/*	This file Test.h is part of scallop.
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

#ifndef SCALLOP_INPUT_TEST_TEST_H_
#define SCALLOP_INPUT_TEST_TEST_H_

#include <string>

namespace scallop {
namespace input {
namespace test {

/**
 * 	Test the input methods of scallop.
 *
 * 	It creates a sample input test file and emulates a call to Setup.
 * 	Then the sample input file is parsed using InputBaseTest.
 */
class Test
{
public:

	/**
	 * 	Run all tests of this suite.
	 */
	void run_test();
private:

	/**
	 * 	Parse a sample input file that is created.
	 */
	void test_input_file_parsing();

	/**
	 * 	Create a sample input file.
	 * @param fileName The file Name of the sample input file.
	 */
	void create_test_input_file(std::string const& fileName);

};

} /* namespace test */
} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_TEST_TEST_H_ */
