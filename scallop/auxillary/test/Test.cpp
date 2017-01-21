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
 *  Created on: Jan 20, 2017
 *      Author: A. Linscheid
 */

#include "Test.h"
#include "scallop/gw_flex/MemoryLayout.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/error_handling/Error.h"
#include <sstream>
#include <complex>
#include <vector>

namespace scallop
{
namespace auxillary
{
namespace test
{

void Test::run_test()
{
	test_NambuSpin_decomposition();
}

void Test::test_NambuSpin_decomposition()
{
	output::TerminalOut msg;

	gw_flex::MemoryLayout mem;
	mem.initialize_layout_2pt_obj(1);
	std::vector< std::complex<double> > nsmat_s_r(16, std::complex<double>(0,0) );
	for ( size_t a1 = 0 ; a1 < 2; ++a1)
		for ( size_t s1 = 0 ; s1 < 2; ++s1)
			for ( size_t a2 = 0 ; a2 < 2; ++a2)
				for ( size_t s2 = 0 ; s2 < 2; ++s2)
				{
					nsmat_s_r[mem.memory_layout_2pt_obj(0,a1,s1,0,a2,s2)] =
							auxillary::BasicFunctions::singlet_Re_channel(a1,s1,a2,s2);
				}
	for ( size_t a1 = 0 ; a1 < 2; ++a1)
		for ( size_t s1 = 0 ; s1 < 2; ++s1)
		{
			std::string s;
			std::stringstream ss(s);
			for ( size_t a2 = 0 ; a2 < 2; ++a2)
				for ( size_t s2 = 0 ; s2 < 2; ++s2)
					ss << nsmat_s_r[mem.memory_layout_2pt_obj(0,a1,s1,0,a2,s2)] << '\t';
			msg << ss.str();
		}

	if ( nsmat_s_r[mem.memory_layout_2pt_obj(0,0,0,0,1,1)] != nsmat_s_r[mem.memory_layout_2pt_obj(0,1,1,0,0,0)] )
		error_handling::Error("Test failed");
}

} /* namespace test */
} /* namespace auxillary */
} /* namespace scallop */
