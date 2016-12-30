/*	This file main.cpp is part of scallop.
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
 *  Created on: Nov 14, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/input/input.h"
#include "scallop/output/output.h"
#include "scallop/error_handling/Warning.h"
#include "scallop/eliashberg/eliashberg.h"
#include "scallop/gw_flex/gw_flex.h"

using namespace scallop;

int main(int argc, char *argv[])
{
	//read the input file //TODO or stdin
	input::Setup setup(argc,argv);

	if ( setup.get_config().get_method().compare("eli") == 0 )
	{
		//solve the Eliashberg equations
		eliashberg::DriverImaginaryAxis<double> eliashEqDriver;

		eliashEqDriver.set_input(setup);

		eliashEqDriver.solve();
	}
	else if ( setup.get_config().get_method().compare("kgw") == 0 )
	{
		//decide the accuracy of the floating point math
		typedef std::complex<double> T;

		gw_flex::Driver<T> flex_driver( setup.get_config() );
		flex_driver.converge();
	}
	else
	{
		error_handling::Error( std::string("unkown task : '")+setup.get_config().get_method()+"'");
	}

	return 0;
};


