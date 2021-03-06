/*	This file Test.cpp is part of scallop.
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

#include "Test.h"

namespace scallop {
namespace eliashberg {
namespace test{

void Test::run_test() {

	//Tests of individual components
	this->run_MatzubaraSplittingVector_test<double>();
	this->run_DriverImaginaryAxis_test<double>();
}

} /* namespace test */
} /* namespace eliashberg */
} /* namespace scallop */
