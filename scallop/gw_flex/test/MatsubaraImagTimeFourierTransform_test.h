/*	This file MatsubaraImagTimeFourierTransform_test.h is part of scallop.
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
 *  Created on: Oct 28, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_TEST_MATSUBARAIMAGTIMEFOURIERTRANSFORM_TEST_H_
#define SCALLOP_GW_FLEX_TEST_MATSUBARAIMAGTIMEFOURIERTRANSFORM_TEST_H_

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class MatsubaraImagTimeFourierTransform_test
{
public:

	void test_free_particle_greensfunction();
};

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/test/MatsubaraImagTimeFourierTransform_test.hpp"
#endif /* SCALLOP_GW_FLEX_TEST_MATSUBARAIMAGTIMEFOURIERTRANSFORM_TEST_H_ */
