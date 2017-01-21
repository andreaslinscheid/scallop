/*	This file GapFileReader_test.hpp is part of scallop.
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
 *  Created on: Jan 18, 2017
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GapFileReader.h"
#include "scallop/error_handling/Error.h"
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

namespace scallop
{
namespace gw_flex
{
namespace test
{

class GapFileReader_test
{
public:
	void test_all();

	void test_dxy_gap();

	void create_test_input_file_dxy_1b(std::string const& filename);

	void create_test_input_file_spm_2b(std::string const& filename);

	void create_test_input_file_pt_1b(std::string const& filename);
};

void GapFileReader_test::test_all()
{
	test_dxy_gap();
}

void GapFileReader_test::test_dxy_gap()
{
	std::string filename = "/tmp/scallop_gap_input.dat";
	create_test_input_file_dxy_1b(filename);

	GapFileReader g;
	g.read_file(filename);
	std::vector<double> ksample1 = { 0.25, 0.25};
	std::vector<double> ksample2 = { 0.00, 0.25};
	std::vector<double> ksample3 = { 0.25, 0.75};
	std::vector<double> ksample4 = { 0.75, 0.75};
	if ( std::abs( std::real(g(ksample1,double(5.5),/* spin-state */0,/* orbital */0)) - (0.5*-0.2*0.54627421529) ) > 0.00001 )
		error_handling::Error("Test failed reading dxy gap");
	if ( std::abs( std::real(g(ksample2,double(5.5),/* spin-state */0,/* orbital */0)) ) > 0.00001 )
		error_handling::Error("Test failed reading dxy gap");
	if ( std::abs( std::real(g(ksample3,double(5.5),/* spin-state */0,/* orbital */0)) + (0.5*-0.2*0.54627421529) ) > 0.00001 )
		error_handling::Error("Test failed reading dxy gap");
	if ( std::abs( std::real(g(ksample4,double(5.5),/* spin-state */0,/* orbital */0)) - (0.5*-0.2*0.54627421529) ) > 0.00001 )
		error_handling::Error("Test failed reading dxy gap");
	if ( std::abs( std::real(g(ksample4,double(5.5),/* spin-state */1,/* orbital */0)) ) > 0.00001 )
		error_handling::Error("Test failed reading dxy gap");

	create_test_input_file_spm_2b(filename);
	g.read_file(filename);
	if ( (std::abs( std::real(g(ksample1,double(100),0,0)) - (0.5*-0.2*0.282094791) ) > 0.00001)
			|| (std::abs( std::real(g(ksample1,double(100),0,1)) - (0.5*0.282094791) ) > 0.00001) )
		error_handling::Error("Test failed reading s pm gap");
	if ( (std::abs( std::real(g(ksample1,double(100),1,0)) ) > 0.00001)
			|| (std::abs( std::real(g(ksample1,double(100),1,1)) ) > 0.00001) )
		error_handling::Error("Test failed reading s pm gap");
}

void GapFileReader_test::create_test_input_file_dxy_1b(std::string const& filename)
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( mpi.ioproc() )
	{
		std::ofstream file( filename.c_str() );
		if ( not file.good() )
			error_handling::Error(std::string("Unable to create file ")+filename);
		file << "s dxy 5.5 -0.2";
	}
}

void GapFileReader_test::create_test_input_file_spm_2b(std::string const& filename)
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( mpi.ioproc() )
	{
		std::ofstream file( filename.c_str() );
		if ( not file.good() )
			error_handling::Error(std::string("Unable to create file ")+filename);
		file << "s s 100.0 -0.2 1.0";
	}
}

void GapFileReader_test::create_test_input_file_pt_1b(std::string const& filename)
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( mpi.ioproc() )
	{
		std::ofstream file( filename.c_str() );
		if ( not file.good() )
			error_handling::Error(std::string("Unable to create file ")+filename);
		file << "tx px 100.0 1.2";
		file << "ty py 100.0 -1.2";
	}
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
