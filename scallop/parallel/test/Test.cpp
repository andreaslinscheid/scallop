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
 *  Created on: Nov 5, 2016
 *      Author: A. Linscheid
 */

#include "scallop/parallel/test/Test.h"
#include "scallop/parallel/test/GridDistribution_test.h"
#include <vector>
#include <complex>

namespace scallop
{
namespace parallel
{
namespace test
{

void Test::run_test()
{
	typedef std::complex<double> T;

	//testing the scallop vector
	auxillary::TemplateTypedefs<T>::scallop_vector aTest(5, T(0) );


	GridDistribution_test<T> grid_test;
	grid_test.simple_grid_3D();
	grid_test.simple_grid_2D();
}

} /* namespace test */
} /* namespace parallel */
} /* namespace scallop */
