/*	This file UnitaryWannierKSBands_test.hpp is part of scallop.
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
 *  Created on: Nov 10, 2016
 *      Author: A. Linscheid
 */

#include "UnitaryWannierKSBands_test.h"
#include "scallop/gw_flex/UnitaryWannierKSBands.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{
template<typename T>
void UnitaryWannierKSBands_test<T>::test_all()
{
	test_initialization();
}

template<typename T>
void UnitaryWannierKSBands_test<T>::test_initialization()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	msg << "Testing UnitaryTransform from KS to Wannier bands";
	UnitaryWannierKSBands<T> unit;
	if ( mpi.get_nproc() < 2 )
	{
		std::vector<size_t> grid = { 1, 1 };
		unit.initialize_identity( grid, 1 );

		T diff = 0;
		for ( size_t m1 = 0 ; m1 < unit.get_nOrb()*4 ; ++m1)
			for ( size_t m2 = 0 ; m2 < unit.get_nOrb()*4 ; ++m2)
			{
				diff += abs(unit(0,m1,m2) - (m1==m2?T(1.0):T(0.0)));
			}
		assert( (diff.real() == 0) && (diff.imag() == 0) );
	}
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
