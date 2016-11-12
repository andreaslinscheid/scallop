/*	This file GridDistribution_test.hpp is part of scallop.
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

#include "scallop/parallel/test/GridDistribution_test.h"
#include "scallop/parallel/GridDistribution.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include <vector>
#include <cstdlib>
#include <cassert>

namespace scallop
{
namespace parallel
{
namespace test
{

template<typename T>
void GridDistribution_test<T>::simple_grid_2D()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	msg << "Testing grid distribution. Currently using "
			<< mpi.get_nproc() << " processors. Do not use more than 3.";
	mpi.barrier();

	std::vector<size_t> testgrid = {3, 5};

	GridDistribution<T> gridD;
	gridD.distribute_grid( testgrid );

	msg << "Testing grid layout to be what we expect (k is row major in x,y,z )";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << " k grid dimension: \n";
			for (auto k:gridD.get_k_grid())
				std::cout << k << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}
	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_to_xyz( ig );
		size_t expect = tuple[0]*gridD.get_k_grid()[1]+tuple[1];
		assert( expect == ig);
	}
	typename auxillary::TemplateTypedefs<T>::scallop_vector data( gridD.get_num_grid_data(), T(0) );

	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
		data[ gridD.k_conseq_local_to_data_conseq(ig) ] = T(tuple[0]+1,tuple[1]+1);
	}

	msg << "Testing data layout per processor, grid is k space 3x5:\n"
			"2D complex data, where the imaginary component is the y coordinate\n"
			<< "Zero padding to meet constant storage space under Fourier transform criteria.";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";
			for (auto d:data)
				std::cout << d << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	msg << "Testing data transposition to R space in column major ordering:\n";
	gridD.grid_data_transposition( true, data, 1 );

	T diff = T(0);
	for ( size_t iR = 0 ; iR < gridD.get_num_R_grid(); ++iR )
	{
		std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( iR );
		diff += std::abs( data[ gridD.R_conseq_local_to_data_conseq(iR) ] - T(tuple[0]+1,tuple[1]+1) );
	}

	mpi.sum(diff);
	msg << "Difference between the grid data after transposition and the return of the object:" << diff << "\n\n";

	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	msg << "Transposing again and comparing to the initial data: ";
	gridD.grid_data_transposition( false, data, 1 );
	msg << "After transposition:\n";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<"):\t"
							<< data[ gridD.k_conseq_local_to_data_conseq(ig) + 0 ]<< "]\n";
			}
			std::cout << std::endl ;
		}
		mpi.barrier();
	}

	msg << "Testing grid layout to be what we expect (R is column major in x,y,z )";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << " R grid dimension: \n";
			for (auto R:gridD.get_R_grid())
				std::cout << R << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}
	for ( size_t ig = 0 ; ig < gridD.get_num_R_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.R_conseq_to_xyz( ig );
		size_t expect = tuple[1]*gridD.get_R_grid()[0]+tuple[0];
		assert( expect == ig);
	}

	msg << "Testing grid in R space 8x8:\n";
	std::vector<size_t> testgrid2 = {8, 8};
	gridD.distribute_grid( testgrid2 );
	data = typename auxillary::TemplateTypedefs<T>::scallop_vector( gridD.get_num_grid_data(), T(0) );
	for ( size_t iR = 0 ; iR < gridD.get_num_R_grid(); ++iR )
	{
		std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( iR );
		data[ gridD.R_conseq_local_to_data_conseq(iR) ] = T(tuple[0]+1,tuple[1]+1);
	}
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";
			for (auto d:data)
				std::cout << d << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}
	gridD.grid_data_transposition( false, data, 1 );
	diff = T(0);
	for ( size_t ik = 0 ; ik < gridD.get_num_R_grid(); ++ik )
	{
		std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ik );
		diff += std::abs( data[ gridD.k_conseq_local_to_data_conseq(ik) ] - T(tuple[0]+1,tuple[1]+1) );
	}
	mpi.sum(diff);
	msg << "Difference between the grid data after transposition and the return of the object:" << diff << "\n\n";

	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	msg << "Testing 10x11 k grid with 2 points of data:\n";
	std::vector<size_t> testgrid3 = {10, 11};
	gridD.distribute_grid( testgrid3 );
	data = typename auxillary::TemplateTypedefs<T>::scallop_vector( gridD.get_num_grid_data()*2, T(0) );
	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
		data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ] = T(tuple[0]+1,tuple[1]+1);
		data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ] = T(tuple[0]+1,tuple[1]+1);
	}

	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<"):\t" << data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ]
				                    << "\t" << data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ] << "]\n";
			}
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	gridD.grid_data_transposition( true, data, 2 );

	msg << "After transposition:\n";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_R_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<"):\t"
							<< data[ gridD.R_conseq_local_to_data_conseq(ig)*2 + 0 ]
					<< "\t" << data[ gridD.R_conseq_local_to_data_conseq(ig)*2 + 1 ] << "]\n";
			}
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	diff = T(0);
	for ( size_t iRy = 0 ; iRy < gridD.get_R_grid()[1]; ++iRy )
		for ( size_t iRx = 0 ; iRx < gridD.get_R_grid()[0]; ++iRx )
		{
			size_t expected_layout = iRy*gridD.get_R_grid()[0]+iRx;
			std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( expected_layout );
			diff += std::abs( data[ gridD.R_conseq_local_to_data_conseq(expected_layout)*2+0 ]
			                        - T(tuple[0]+1,tuple[1]+1) );
			diff += std::abs( data[ gridD.R_conseq_local_to_data_conseq(expected_layout)*2+1 ]
			                        - T(tuple[0]+1,tuple[1]+1) );
		}
	mpi.sum(diff);
	msg << "Computing difference to the expected result: " << diff << "\n\n";

	mpi.barrier();
	assert( ( diff.real() < 0.00001 )&&( diff.imag() < 0.00001 ) );

	msg << "Transposing again and comparing to the initial data: ";
	gridD.grid_data_transposition( false, data, 2 );
	msg << "After transposition:\n";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<"):\t"
							<< data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ]
					<< "\t" << data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ] << "]\n";
			}
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
		diff += std::abs( std::real(data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ])
						- std::real(T(tuple[0]+1,tuple[1]+1)));
		diff += std::abs( std::real(data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ])
						- std::real(T(tuple[0]+1,tuple[1]+1)));

		diff += std::abs( std::imag(data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ])
						- std::imag(T(tuple[0]+1,tuple[1]+1)));
		diff += std::abs( std::imag(data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ])
						- std::imag(T(tuple[0]+1,tuple[1]+1)));
	}

	mpi.sum(diff);
	msg << "Computing difference to the expected result: " << diff << "\n\n";

	mpi.barrier();
	assert( ( diff.real() < 0.00001 )&&( diff.imag() < 0.00001 ) );
}

template<typename T>
void GridDistribution_test<T>::simple_grid_3D()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	msg << "Testing grid distribution. Currently using "
			<< mpi.get_nproc() << " processors. Do not use more than 3.";
	mpi.barrier();
	std::vector<size_t> testgrid = {5,4,3};
	GridDistribution<T> gridD;
	gridD.distribute_grid( testgrid );

	msg << "Testing 3D grid layout to be what we expect (k is row major in x,y,z )";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << " k grid dimension: \n";
			for (auto k:gridD.get_k_grid())
				std::cout << k << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}
	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_to_xyz( ig );
		size_t expect = (tuple[0]*gridD.get_k_grid()[1]+tuple[1])*gridD.get_k_grid()[2]+tuple[2];
		assert( expect == ig);
	}
	for ( size_t ig = 0 ; ig < gridD.get_num_R_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.R_conseq_to_xyz( ig );
		size_t expect = (tuple[2]*gridD.get_R_grid()[1]+tuple[1])*gridD.get_R_grid()[0]+tuple[0];
		assert( expect == ig);
	}


	msg << "Testing 5x4x3 k grid with 2 points of data:\n";
	auto data = typename auxillary::TemplateTypedefs<T>::scallop_vector( gridD.get_num_grid_data()*2, T(0) );
	for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
	{
		std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
		data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ] = T(tuple[0]+1,tuple[1]+1);
		data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ] = T(tuple[0]+1,tuple[1]+1);
	}

	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_k_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.k_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<",iz="<<tuple[2]<<"):\t"
						<< data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 0 ]
						<< "\t" << data[ gridD.k_conseq_local_to_data_conseq(ig)*2 + 1 ] << "]\n";
			}
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	gridD.grid_data_transposition( true, data, 2 );

	msg << "After transposition:\n";
	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << ": \n";

			for ( size_t ig = 0 ; ig < gridD.get_num_R_grid(); ++ig )
			{
				std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( ig );
				std::cout << "[(ix="<< tuple[0] << ",iy="<<tuple[1]<<",iz="<<tuple[2]<<"):\t"
							<< data[ gridD.R_conseq_local_to_data_conseq(ig)*2 + 0 ]
					<< "\t" << data[ gridD.R_conseq_local_to_data_conseq(ig)*2 + 1 ] << "]\n";
			}
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	T diff = T(0);
	for ( size_t iRz = 0 ; iRz < gridD.get_R_grid()[2]; ++iRz )
		for ( size_t iRy = 0 ; iRy < gridD.get_R_grid()[1]; ++iRy )
			for ( size_t iRx = 0 ; iRx < gridD.get_R_grid()[0]; ++iRx )
			{
				size_t expected_layout = (iRz*gridD.get_R_grid()[1]+iRy)*gridD.get_R_grid()[0]+iRx;
				std::vector<size_t> tuple = gridD.R_conseq_local_to_xyz_total( expected_layout );
				diff += std::abs( data[ gridD.R_conseq_local_to_data_conseq(expected_layout)*2+0 ]
										- T(tuple[0]+1,tuple[1]+1) );
				diff += std::abs( data[ gridD.R_conseq_local_to_data_conseq(expected_layout)*2+1 ]
										- T(tuple[0]+1,tuple[1]+1) );
			}
	mpi.sum(diff);
	msg << "Computing difference to the expected result: " << diff << "\n\n";
}

} /* namespace test */
} /* namespace parallel */
} /* namespace scallop */
