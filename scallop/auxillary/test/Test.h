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
 *  Created on: Jan 20, 2017
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_AUXILLARY_SRC_TEST_H_
#define SCALLOP_AUXILLARY_SRC_TEST_H_

namespace scallop
{
namespace auxillary
{
namespace test
{

class Test
{
public:

	void run_test();

private:

	void test_NambuSpin_decomposition();
};

} /* namespace test */
} /* namespace auxillary */
} /* namespace scallop */

#endif /* SCALLOP_AUXILLARY_SRC_TEST_H_ */
