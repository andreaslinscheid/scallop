/*	This file GridDistribution_test.h is part of scallop.
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
 *  Created on: Nov 5, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_PARALLEL_TEST_GRIDDISTRIBUTION_TEST_H_
#define SCALLOP_PARALLEL_TEST_GRIDDISTRIBUTION_TEST_H_

namespace scallop
{
namespace parallel
{
namespace test
{

template<typename T>
class GridDistribution_test
{
public:

	void simple_grid_2D();

	void simple_grid_3D();
};

} /* namespace test */
} /* namespace parallel */
} /* namespace scallop */

#include "scallop/parallel/test/GridDistribution_test.hpp"
#endif /* SCALLOP_PARALLEL_TEST_GRIDDISTRIBUTION_TEST_H_ */
