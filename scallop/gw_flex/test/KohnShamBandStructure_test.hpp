/*	This file KohnShamBandStructure_test.hpp is part of scallop.
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
 *  Created on: Nov 25, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/KohnShamBandStructure.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class KohnShamBandStructure_test
{
public:

	void test_all();

private:

	void test_2_bnd_cos();
};

template<typename T>
void KohnShamBandStructure_test<T>::test_all()
{
	this->test_2_bnd_cos();
}

template<typename T>
void KohnShamBandStructure_test<T>::test_2_bnd_cos()
{
	KohnShamBandStructure<T> e;

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	const std::string filename = "/tmp/2bndcos.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream testFile( filename.c_str() );
		if ( ! testFile.good() )
			error_handling::Error( std::string()+"Unable to create "+filename+" for testing the input of a 2 band model");
		testFile << "model\nTwoBandCosine\nt1=50.9\nt2=100\nEe=-20\nEh=10.0\n";
		testFile.close();
	}
	e.initialize_from_file( {64, 64}, filename );

	auto grid = e.get_spaceGrid_proc().get_grid();
	for (size_t ik = 0 ; ik < e.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = e.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = e.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		size_t mid = grid[1]*grid[0]/2+grid[0]/2;
		if ( ictotal == mid)
		{
			std::cout << "\tValue at Q=(pi,pi): band 1: " << e(ik,0,0,0)
						<< "eV; band 2: " << e(ik,1,0,0) << "eV"  << std::endl;
			assert( (std::abs( e(ik,1,0,0) + 20.0) < 0.0000001) );
		}
		if ( ictotal == 0 )
		{
			std::cout << "\tValue at Q=(0,0): band 1: " << e(0,0,0,0)
						<< "eV; band 2: " << e(0,1,0,0) <<"eV" << std::endl;
			assert( (std::abs( e(0,0,0,0) - 10.0) < 0.0000001) );
		}
	}
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
