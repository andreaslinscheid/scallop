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
 *  Created on: Nov 25, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_TEST_H_
#define SCALLOP_ELIASHBERG_TEST_H_

#include "scallop/eliashberg/eliashberg.h"

namespace scallop {
namespace eliashberg {
namespace test{

class Test {
public:
	void run_test();
private:

	template<typename T>
	void run_MatzubaraSplittingVector_test();

	template<typename T>
	void run_DriverImaginaryAxis_test();
};

} /* namespace test */
} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/test/Test.hpp"
#endif /* SCALLOP_ELIASHBERG_TEST_H_ */
