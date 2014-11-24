/*	This file Setup.cpp is part of scallop.
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

#include "scallop/input/Setup.h"
#include "scallop/error_handling/Error.h"
#include <fstream>
#include <vector>
#include <iostream>

namespace scallop {
namespace input {

Setup::Setup(int argc, char *argv[]) {
	//Print header

	//If we are in debugging mode, generate a stack trace in the case of a segmentation fault
	error_handling::Error::generate_stacktrace_on_segfault();

	if ( argc < 1 ) {
		error_handling::Error( "\tNo options passed, exiting ... !", 1 );
	}else if ( argc == 2 ) {
		//fill the internal key/value map with content.
		_inputFile.read_input_file(argv[1]);
	} else {
		error_handling::Error( "\tCannot figure out how to get input, exiting ... !", 1 );
	}
}


InputFile const& Setup::get_parsed_input_file() const {
	return _inputFile;
}

std::vector<std::string> Setup::get_list_unread_input_parameters() const {
	return _inputFile.get_list_unread_input_parameters();
}

} /* namespace input */
} /* namespace scallop */
