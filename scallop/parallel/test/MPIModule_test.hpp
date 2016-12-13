/*	This file MPIModule_test.hpp is part of scallop.
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
 *  Created on: Nov 30, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_PARALLEL_MPIMODULE_TEST_H_
#define SCALLOP_PARALLEL_MPIMODULE_TEST_H_

#include "scallop/parallel/MPIModule.h"
#include <mpi.h>
#include <complex>
#include <vector>

namespace scallop
{
namespace parallel
{
namespace test
{

template<typename T>
class MPIModule_test
{
public:
	void test_all();

	void test_alltoallv();
};

template<typename T>
void MPIModule_test<T>::test_all()
{
	this->test_alltoallv();
}

template<typename T>
void MPIModule_test<T>::test_alltoallv()
{
	MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( mpi.get_nproc() == 2 )
	{
		size_t me = mpi.get_mpi_me();

		//On proc 0 we put in 4 elements and on proc 1, 2
		std::vector< T > data;
		if ( me == 0 )
			data = { T(1,me), T(2,me), T(3,me), T(4,me)};
		if ( me == 1 )
			data = { T(1,me), T(2,me)};

		std::vector<size_t> sendc( mpi.get_nproc() );
		//The first 2 stay on this proc, the second
		if ( me == 0 )
			sendc = {2,2};
		//The 1 are send to the first this proc
		if ( me == 1 )
			sendc = {1,1};

		std::vector< T > rec;
		mpi.all_to_allv( data, rec, sendc );

		//On proc 0 we have copied the first two elements and
		//	received the 3 elements from proc 1
		std::vector< T > expect;
		if ( me == 0 )
			expect = { T(1,0), T(2,0), T(1,1) };
		if ( me == 1 )
			expect = { T(3,0), T(4,0), T(2,1) };

		if ( ! (expect == rec) )
			error_handling::Error( std::string("Check failed on proc ")+std::to_string(me));
	}
}

} /* namespace test */
} /* namespace parallel */
} /* namespace scallop */

#endif /* SCALLOP_PARALLEL_MPIMODULE_TEST_H_ */
