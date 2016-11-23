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
 *  Created on: Oct 28, 2016
 *      Author: alinsch
 */

#ifndef SCALLOP_GW_FLEX_TEST_TEST_H_
#define SCALLOP_GW_FLEX_TEST_TEST_H_

namespace scallop {
namespace gw_flex {
namespace test {

class Test
{
public:
	void run_test();

private:

	void test_time_freq_Fourier_transform();

	void test_GreensFunctionOrbital();

	void test_FFTBase();

	void test_UnitaryWannierKSBands();

	void test_Susceptibility();

	void test_phonon_gf();

	void test_InteractionMatrix();
};

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */

#endif /* SCALLOP_GW_FLEX_TEST_TEST_H_ */
