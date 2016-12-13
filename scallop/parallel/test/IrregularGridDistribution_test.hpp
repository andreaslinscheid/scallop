/*	This file IrregularGridDistribution_test.hpp is part of scallop.
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
 *  Created on: Nov 28, 2016
 *      Author: A. Linscheid
 */

#include "scallop/parallel/IrregularGridDistribution.h"
#include "scallop/parallel/GridDistribution.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/output/TerminalOut.h"

namespace scallop
{
namespace parallel
{
namespace test
{

template<typename T>
class IrregularGridDistribution_test
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	void test_all();

private:

	void test_construction();

	void test_data_sync();
};

template<typename T>
void IrregularGridDistribution_test<T>::test_all()
{
	this->test_construction();

	this->test_data_sync();
}

template<typename T>
void IrregularGridDistribution_test<T>::test_construction()
{
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	if ( mpi.get_nproc() == 1 )
	{
		msg << "\nConstructing a 5x5 regular grid and connecting the"
				"point (0.01,0.1)";
		typename auxillary::TemplateTypedefs<bT>::scallop_vector testPts = {0.01,0.1};
		GridDistribution<T> gdist;
		gdist.distribute_grid( {5,5});

		IrregularGridDistribution<T> someIrregularPts;
		someIrregularPts.distribute_pts( /*inKSpace=*/true,	testPts,gdist);

		//the pt is in (0.0,0.0) (0.2,0.2)
		std::vector<size_t> listOfIndices;
		someIrregularPts.get_list_required_regular_grid_indices(listOfIndices);

		std::vector<size_t> expected = {0,1,5,6};
		if ( listOfIndices != expected )
			error_handling::Error("Check failed.");

		V dataFullGrid(gdist.get_num_k_grid(),T(0));
		for ( size_t ik = 0 ; ik < gdist.get_num_k_grid(); ++ik)
		{
			auto tuple = gdist.k_conseq_local_to_xyz_total(ik);
			dataFullGrid[ik] = T(tuple[0]+1,tuple[1]+1);
		}

		V dataSelection;
		for (auto i : listOfIndices)
			dataSelection.push_back( dataFullGrid[i] );

		V interpol;
		someIrregularPts.linear_interpolate_data(dataSelection,interpol,1);
		if ( std::abs(interpol[0] - T(1.05,1.5) ) > 0.000001 )
			error_handling::Error("Check failed.");

		msg << "\t connecting the points (0,0.0),(0.1,0.7) and (-0.1,1.0) with this grid";

		testPts = {0,0.0,0.1,0.7,-0.1,1.0};

		someIrregularPts.distribute_pts( /*inKSpace=*/true,	testPts,gdist);
		someIrregularPts.get_list_required_regular_grid_indices(listOfIndices);
		expected = {0,1,3,4,5,6,8,9,20,21};
		if ( listOfIndices != expected )
			error_handling::Error("Check failed.");

		dataSelection.clear();
		for (auto i : listOfIndices)
			dataSelection.push_back( dataFullGrid[i] );

		interpol.clear();
		someIrregularPts.linear_interpolate_data(dataSelection,interpol,1);
		if ( std::abs(interpol[0] - T(1.0,1.0) ) > 0.000001 )
			error_handling::Error("Check failed.");
		if ( std::abs(interpol[1] - T(1.5,4.5) ) > 0.000001 )
			error_handling::Error("Check failed.");
		if ( std::abs(interpol[2] - T(3.0,1.0) ) > 0.000001 )
			error_handling::Error("Check failed.");
	}
	else if ( mpi.get_nproc() == 2 )
	{
		msg << "\nConstructing a 5x5 regular grid and connecting"
				"the points (0,0.0),(0.1,0.7) and (-0.1,1.0) and with 2 processors";
		typename auxillary::TemplateTypedefs<bT>::scallop_vector testPts
						= {0,0.0,0.1,0.7,-0.1,1.0};
		GridDistribution<T> gdist;
		gdist.distribute_grid( {5,5});

		IrregularGridDistribution<T> someIrregularPts;
		someIrregularPts.distribute_pts( /*inKSpace=*/true,	testPts,gdist);

		std::vector<size_t> listOfIndices;
		someIrregularPts.get_list_required_regular_grid_indices(listOfIndices);

		//kx is parallelized, proc 0 handles 0,1,2 and proc 1 handles 3,4
		if ( mpi.get_mpi_me() == 0 )
		{
			std::vector<size_t> expected = {0,1,3,4,5,6,8,9};
			if ( listOfIndices != expected )
			{
				for (auto a : listOfIndices )
					msg << '\t'<< a ;
				error_handling::Error("Check failed on proc 0.");
			}
		}
		if ( mpi.get_mpi_me() == 1 )
		{
			std::vector<size_t> expected = {0,1,20,21};
			if ( listOfIndices != expected )
				error_handling::Error("Check failed on proc 1.");
		}
	}
}

template<typename T>
void IrregularGridDistribution_test<T>::test_data_sync()
{
	MPIModule const& mpi = parallel::MPIModule::get_instance();
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	output::TerminalOut msg;

	msg << "\nTesting data synchronization of an irregular with a regular grid";
	msg << "\tRegular grid of 5x5 and 3 points";
	GridDistribution<T> gdist;
	gdist.distribute_grid( {5,5} );
	typename auxillary::TemplateTypedefs<bT>::scallop_vector testPts = {0,0.0,0.1,0.7,-0.1,1.0};

	V dataProcGrid(gdist.get_num_k_grid(),T(0));
	for ( size_t ik = 0 ; ik < gdist.get_num_k_grid(); ++ik)
	{
		auto tuple = gdist.k_conseq_local_to_xyz_total(ik);
		dataProcGrid[ik] = T(tuple[0]+1,tuple[1]+1);
	}

	bool inKSpace = true;
	IrregularGridDistribution<T> someIrregularPts;
	someIrregularPts.distribute_pts( inKSpace, testPts, gdist);

	if (mpi.get_mpi_me() == 0)
	{
		std::cout << "\tProcessor 0 has the grid points:" << std::endl;
		for ( size_t i = 0 ; i < someIrregularPts.get_n_pts(); ++i)
		{
			for ( auto ki : someIrregularPts.get_vector(i))
				std::cout << '\t'<< ki ;
			std::cout << std::endl;
		}
	}
	mpi.barrier();
	if (mpi.get_mpi_me() == 1)
	{
		std::cout << "\tProcessor 1 has the grid points:" << std::endl;
		for ( size_t i = 0 ; i < someIrregularPts.get_n_pts(); ++i)
		{
			for ( auto ki : someIrregularPts.get_vector(i))
				std::cout <<'\t' << ki;
			std::cout << std::endl;
		}
	}
	mpi.barrier();

	msg << "\tTesting multi-processor interpolation on irregular grid points:";
	V dataForThisPts;
	someIrregularPts.proc_sync_data( inKSpace, dataProcGrid,  dataForThisPts, 1 );

	V interpol;
	someIrregularPts.linear_interpolate_data(dataForThisPts,interpol,1);

	for ( size_t i = 0 ; i < someIrregularPts.get_n_pts() ; ++i)
	{
		size_t global = someIrregularPts.proc_local_to_total_index( i );
		if ( global == 0 )
			if ( std::abs(interpol[i] - T(1.0,1.0) ) > 0.000001 )
				error_handling::Error("Check 1 failed.");
		if ( global == 1 )
			if ( std::abs(interpol[i] - T(1.5,4.5) ) > 0.000001 )
				error_handling::Error("Check 2 failed.");
		if ( global == 2 )
			if ( std::abs(interpol[i] - T(3.0,1.0) ) > 0.000001 )
				error_handling::Error("Check 3 failed.");
	}
	msg << "\tcomplete.";


	msg << "\nTesting data synchronization of an irregular with a regular grid";
	msg << "\tRegular grid of 4x4x2 and 3 points, with 2 data values";
	gdist.distribute_grid( {4,4,2} );
	dataProcGrid = V(gdist.get_num_k_grid()*2,T(0));
	for ( size_t ik = 0 ; ik < gdist.get_num_k_grid(); ++ik)
	{
		size_t totalIndex = gdist.k_xyz_to_conseq(gdist.k_conseq_local_to_xyz_total(ik));
		dataProcGrid[ik*2] = T(totalIndex,std::pow(2,mpi.get_mpi_me()+1));
		dataProcGrid[ik*2+1] = T(totalIndex,std::pow(2,mpi.get_mpi_me()+1)+1);
	}
	testPts = {0.0, 0.0, 0.0,
			   0.1, 0.7, 0.0,
			  -0.1, 1.0, 0.0};
	someIrregularPts.distribute_pts( inKSpace, testPts, gdist);
	dataForThisPts = V();
	someIrregularPts.proc_sync_data( inKSpace, dataProcGrid,  dataForThisPts, 2 );
	if (mpi.get_mpi_me() == 0)
	{
		std::cout << "\tProcessor 0 has the data :" << std::endl;
		for ( auto d : dataForThisPts)
			std::cout << '\t'<< d ;
		std::cout << std::endl;
	}
	mpi.barrier();
	if (mpi.get_mpi_me() == 1)
	{
		std::cout << "\tProcessor 1 has the data :" << std::endl;
		for ( auto d : dataForThisPts)
			std::cout << '\t'<< d ;
	}
	mpi.barrier();
}

} /* namespace test */
} /* namespace parallel */
} /* namespace scallop */
