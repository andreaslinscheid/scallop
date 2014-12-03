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
#include "scallop/error_handling/Warning.h"
#include "scallop/eliashberg/eliashberg.h"

using namespace scallop;

int main(int argc, char *argv[]) {

	//read the input file //TODO or stdin
	input::Setup setup(argc,argv);

	//Set basic configuration options
	input::Configuration config(setup.get_parsed_input_file());

	if ( config.get_method().compare("Eliash") == 0 ){

		//solve the Eliashberg equations
		eliashberg::DriverImaginaryAxis<double> eliashEqDriver;
		eliashEqDriver.set_input(setup);
		do {
			eliashEqDriver.iterate();
		} while ( not eliashEqDriver.converged() );
	}

	//Flag a warning for all keys that have not been used, as this indicates unintended input
	std::vector<std::string> const unusedKeys = setup.get_list_unread_input_parameters();
	for ( auto & key : unusedKeys )
		error_handling::Warning("Key "+ key +" defined on input has never been used!") ;

	return 0;
};


