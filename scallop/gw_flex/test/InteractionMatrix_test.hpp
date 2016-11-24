/*	This file InteractionMatrix_test.hpp is part of scallop.
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
 *  Created on: Nov 22, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_INTERACTIONMATRIX_TEST_H_
#define SCALLOP_GW_FLEX_INTERACTIONMATRIX_TEST_H_

#include "scallop/gw_flex/InteractionMatrix.h"
#include <fstream>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class InteractionMatrix_test
{
public:
	void test_all();
	void test_init_file();
private:
};

template<typename T>
void InteractionMatrix_test<T>::test_all()
{
	this->test_init_file();
}


template<typename T>
void InteractionMatrix_test<T>::test_init_file()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	const std::string fnameSpin = "/tmp/test_input_s.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream testFile( fnameSpin.c_str() );
		//Testing InteractionMatrix input from file:
		testFile << "1 4" << '\n'; //One orbital 4 channels
		testFile << "0\t0\t0\t0\t0\t0\t1.5\t0.0" << '\n';
		testFile << "1\t1\t0\t0\t0\t0\t3.0\t0.0" << '\n';
		testFile << "2\t2\t0\t0\t0\t0\t4.5\t0.0" << '\n';
		testFile << "3\t3\t0\t0\t0\t0\t6.0\t0.0" << '\n';
		//			^j ^jp l1 l2 l3 l4 Re(I) Im(I)
		testFile.close();
	}

	const std::string fnameCharge = "/tmp/test_input_c.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream testFile( fnameCharge.c_str() );
		//Testing InteractionMatrix input from file:
		testFile << "1 1" << '\n'; //One orbital one channel
		testFile << "0\t0\t0\t0\t0\t0\t9.0\t0.0" << '\n';
		//			^j ^jp l1 l2 l3 l4 Re(I) Im(I)
		testFile.close();
	}

	InteractionMatrix<T> interact;
	interact.init_file( fnameSpin );

	assert( interact.get_nOrb() == 1 );
	size_t nO = interact.get_nOrb();

#ifndef NDEBUG
	for (size_t j = 0 ; j < 4; j++)
		for (size_t jp = 0 ; jp < 4; jp++)
			for (size_t l1 = 0 ; l1 < nO; l1++)
				for (size_t l2 = 0 ; l2 < nO; l2++)
					for (size_t l3 = 0 ; l3 < nO; l3++)
						for (size_t l4 = 0 ; l4 < nO; l4++)
						{
							T val = j==jp? 1.5*T(j+1,0) : T(0);
							assert( interact(j,jp,l1,l2,l3,l4) == val);
						}
#endif
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
#endif /* SCALLOP_GW_FLEX_INTERACTIONMATRIX_TEST_H_ */
